"""
Classes and functions to facilitate interaction with a PBS queue.

Copyright (C) 2014 University of Southern California, Jie Ren,
Emad Bahrami Samani, Philip J. Uren, Andrew D. Smith

Authors: Jie Ren, Emad Bahrami Samani, Philip J. Uren, Andrew D. Smith

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.

"""

# standard python imports
import time
import unittest
import sys
import StringIO
import getpass
import subprocess32

# for mocking during unit testing
import mock

# pyokit imports
from pyokit.common.pyokitError import PyokitError


###############################################################################
#                                CONSTANTS                                    #
###############################################################################

# Default max vals for queue resource usage
DEF_Q_JOB_MAX = 300
DEF_Q_NODES_MAX = 300
DEF_Q_PROCESSOR_MAX = 300
DEF_Q_MEM_MAX_GB = 1000

# controlling status update frequency
DEFAULT_MIN_THINK_TIME_SECONDS = 300


###############################################################################
#                            EXCEPTION CLASSES                                #
###############################################################################

class PBSDuplicateQueueError(PyokitError):
  pass


class PBSNoSuchQueueError(PyokitError):
  pass


class PBSTrackerError(Exception):

  """Raised for jobs that exceed maximum constraints."""

  def __init__(self, msg):
    """Constructor for PBS constraint violation errors."""
    self.value = msg

  def __str__(self):
    """:return: string representation of this error."""
    return repr(self.value)


class PBSConstraintError(Exception):

  """Raised for jobs that exceed maximum constraints."""

  def __init__(self, msg):
    """Constructor for PBS constraint violation errors."""
    self.value = msg

  def __str__(self):
    """:return: string representation of this error."""
    return repr(self.value)


###############################################################################
#                        CLASSES REPRESENTING PBS JOBS                        #
###############################################################################

class PBSJob(object):

  """Represents a PBS job that has been submitted via this interface."""

  def __init__(self, job_name, walltime_requested, nodes_requested,
               processors_per_node_requested, mem_requested, vmem_requested,
               pmem_requested, command, job_id=None):
    """Constructor for PBS job; see class docstring for param details."""
    self.name = job_name
    self.walltime_requested = walltime_requested
    self.nodes_requested = nodes_requested
    self.processors_per_node_requested = processors_per_node_requested
    self.mem_requested = mem_requested
    self.vmem_requested = vmem_requested
    self.pmem_requested = pmem_requested
    self.command = command
    self.id = job_id

  @property
  def total_processors_requested(self):
    """How many processors, across all nodes, were requested by the job."""
    return self.nodes_requested * self.processors_per_node_requested


class _PBSRunningJob(object):

  """Represents a PBS job that is currently executing."""

  def __init__(self, name, user, time_used, status, queue_name, job_id):
    self.name = name
    self.user = user
    self.time_used = time_used
    self.status = status
    self.queue_name = queue_name
    self.id = job_id


###############################################################################
#        FUNCTIONS CLASSES FOR TRACKING RESOURCE USAGE IN PBS QUEUES          #
###############################################################################

class PBSQueueResourceStatus(object):

  """
  Representation of resources usage for a queue.

  :param job_count:        number of running jobs
  :param nodes_count:      number of nodes being used.
  :param processors_count: number of processors being used.
  :param mem_count:        total amount of memory used
  """

  def __init__(self, job_count=DEF_Q_JOB_MAX, nodes_count=DEF_Q_NODES_MAX,
               processors_count=DEF_Q_PROCESSOR_MAX,
               mem_count_gb=DEF_Q_MEM_MAX_GB):
    """
    Constructor for PBSQueueResourceConstraint.

    See class docstring for parameter definitions. Default values are max vals,
    so calling constructor with no arguments gives a default maxmimum
    constraint on a queue.
    """
    self.job_count = job_count
    self.nodes_count = nodes_count
    self.processors_count = processors_count
    self.mem_count_gb = mem_count_gb


###############################################################################
#                      CLASSES FOR MANAGING PBS QUEUES                        #
###############################################################################

class _PBSTracker(object):

  def __init__(self):
    self.queues = {}

  def __get_queue(self, name):
    if name not in self.queues:
      raise PBSNoSuchQueueError("No such PBS queue: " + name)
    else:
      return self.queues[name]

  def register_queue(self, name, constraints, allow_qstat=False,
                     loud_qstat=True):
    if name in self.queues:
      raise PBSDuplicateQueueError("Queue with name " + name +
                                   " already exists")
    self.queues[name] = _QueueTracker(name, constraints, allow_qstat,
                                      loud_qstat)

  def push_job(self, job, queue_name, verbose=False):
    """Create a new job and push it to the specified queue."""
    q = self.__get_queue(queue_name)
    q.submit(job, verbose=verbose)

  def is_complete(self, jid, queue_name=None):
    if queue_name is not None:
      q = self.__get_queue(queue_name)
      q.is_job_complete(jid)
    else:
      raise PyokitError("oops, not implemented yet")


class _QueueTracker(object):

  """
  Abstract representation of a PBS queue.

  :param q_name:       name of the queue; used when submitting jobs
  :param q_constraint: resource constraints for this queue. If exceeded,
                       submission of new jobs will block until resources
                       become free (if use of qstat is permitted), or raise
                       an exception (if use of qstat is not permitted).
  :param allow_qstat:  if True, allow use of qstat to monitor job completion in
                       this queue. Otherwise, no monitoring is done and all
                       jobs are assumed to execute indefinately.
  :param loud_qstat:   if True, any calls to qstat as a sub-process are
                       notified to the user via stderr
  """

  def __init__(self, q_name, q_constraint=None, allow_qstat=False,
               loud_qstat=True):
    """Init a queue tracker. See class docstring for param descriptions."""
    self._running_jobs = {}
    self._submitted_jobs = {}
    self._queue_name = q_name
    self._allow_qstat = allow_qstat
    self._last_think = None
    self._loud_qstat = loud_qstat

    # for tracking how much queue resources we're using
    self._username = getpass.getuser().strip()
    self._min_think_time_seconds = DEFAULT_MIN_THINK_TIME_SECONDS
    self._resource_constraints = PBSQueueResourceStatus()
    self._resource_status = None
    self.__think()

  def is_running(self, jid):
    self.__think()
    return jid in self._running_jobs

  def nodes_in_use(self):
    """:return: number of nodes currently used by jobs running or queued."""
    self.__think()
    return self._resource_status.nodes_count

  def processors_in_use(self):
    """:return: number of processors used by jobs running or queued."""
    self.__think()
    return self._resource_status.processors_count

  def jobs_running_or_queued(self):
    """:return: a list of job IDs for jobs that are currently running."""
    self.__think()
    return self._running_jobs.keys()

  def num_jobs_running_or_queued(self):
    """:return: number of jobs currently running or waiting to run."""
    self.__think()
    return self._resource_status.job_count

  def __think(self):
    """
    Determine which jobs are currently running and update resource usage stats.

    If qstat is not allowed, we just assume all jobs run indefinately, and that
    no jobs were already running for the user on this queue.
    """
    if not self._allow_qstat:
      if self._resource_status is None:
        # have to assume nothing has been submitted to the queue yet..
        self._resource_status =\
            PBSQueueResourceStatus(job_count=0, nodes_count=0,
                                   processors_count=0)
    else:
      current_time = time.time()
      if (self._last_think is None or
         (current_time - self._last_think) > self._min_think_time_seconds):
        if self._loud_qstat:
          sys.stderr.write("querying queue status... ")
        self._last_think = current_time

        # Open a pipe to the qstat command.
        output = subprocess32.check_output(["qstat", self._queue_name],
                                           universal_newlines=True)

        self._running_jobs = {}
        for l in output.split("\n"):
          l = l.strip()
          if l == "":
            continue  # blank line, skip it.

          f = l.split()[0].strip().lower()
          if f.lower() == "job" or (len(set(f)) == 1 and (f[0] == "-")):
            continue  # header line, or formatting line, skip it.

          jid, nm, user, time_used, status, queue = l.split()

          if queue != self._queue_name or user != self._username:
            continue  # wrong queue or user; skip it.

          job = _PBSRunningJob(nm, user, time_used, status, queue, jid)
          assert(jid not in self._running_jobs)
          self._running_jobs[jid] = job

        # update rescource usage stats
        node_c = 0
        proc_c = 0
        for jid in self._running_jobs:
          node_c += (self._submitted_jobs[jid].nodes_requested
                     if jid in self._submitted_jobs else 1)
          proc_c += (self._submitted_jobs[jid].processors_requested
                     if jid in self._submitted_jobs else 1)
        self._resource_status =\
            PBSQueueResourceStatus(len(self._running_jobs), node_c, proc_c)

  def __has_capacity(self, job):
    if (self.num_jobs_running_or_queued() + 1 >
       self._resource_constraints.job_count):
      return False
    if (self.nodes_in_use() + job.nodes_requested >
       self._resource_constraints.nodes_count):
      return False
    if (self.processors_in_use() + job.total_processors_requested >
       self._resource_constraints.processors_count):
      return False
    # if self.running_mem + job.mem > ribocopSettings.max_pbs_mem:
    #  return False
    # if self.running_vmem + job.vmem > ribocopSettings.max_pbs_vmem:
    #  return False
    # if self.running_pmem + job.pmem > ribocopSettings.max_pbs_pmem:
    #  return False
    return True

  def submit(self, job, verbose=False):
    """
    Submit a job to this queue.

    Submission will block if there isn't enough capacity in this queue to
    accept the job, or another job was submitted to this queue too recently.
    Once capacity is freed up, and/or the cool-down time elapses, the job will
    be submitted.

    :param job: the job ID will be assigned to job.id
    """
    # if we don't have capcity and we're not allow to check status, give up
    if not self.__has_capacity(job) and not self.allow_qstat:
      raise PBSConstraintError("Not enough resources")

    # check whether we have enough spare capacity to submit this job. block if
    # not and wait for queue to free up a bit.
    while not self.__has_capacity(job):
      self.__think()

    # prepare job string for qsub
    job_string = """#!/bin/bash
    #PBS -q %s
    #PBS -N %s
    #PBS -l walltime=%s
    #PBS -l nodes=%i:ppn=%i
    #PBS -l mem=%iM
    #PBS -l vmem=%iM
    #PBS -l pmem=%iM
    %s""" % (self._queue_name, job.name, job.walltime_requested,
             job.nodes_requested, job.processors_per_node_requested,
             job.mem_requested, job.vmem_requested, job.pmem_requested,
             job.command)

    # Open a pipe to the qsub command.
    p = subprocess32.Popen("qsub", universal_newlines=True,
                           stdin=subprocess32.PIPE,
                           stdout=subprocess32.PIPE, close_fds=True)
    try:
      outs, errs = p.communicate(input=job_string, timeout=15)
      job.id = outs.read().strip()
      assert(job.id not in self._running_jobs)
      self._running_jobs[job.id] = job
      if verbose:
        sys.stderr.write("qsub accepted job with JID " + str(job.id))
    except subprocess32.TimeoutExpired:
      p.kill()
      outs, errs = p.communicate()


###############################################################################
#                         UNIT TESTS FOR THIS MODULE                          #
###############################################################################

def build_popen_mock_side_effects(results):
  """
  Build a mock side effects function for popen.

  :param results: a dictionary where keys are process names and values are
                  lists of results; the nth item in a list being the value to
                  provide for stdout on the nth call of that process.
  """
  pname_call_count = {}

  def popen_mock_side_effect(*args, **kwargs):
    class PopenDummy:
      def __init__(self):
        self.p_name = args[0]
        if self.p_name not in pname_call_count:
          pname_call_count[self.p_name] = 0
        self.results = results

      def communicate(self, *args, **kwargs):
        calls = pname_call_count[self.p_name]
        if self.p_name in self.results:
          if calls < len(results[self.p_name]):
            o_stng = self.results[self.p_name][calls]
            stdout_strm = StringIO.StringIO(o_stng)
            stderr_strm = StringIO.StringIO()
            pname_call_count[self.p_name] += 1
            return stdout_strm, stderr_strm
          else:
            raise ValueError("expected call to " + self.p_name)
        else:
          raise ValueError("unexpected call to popen")
    return PopenDummy()
  return popen_mock_side_effect


class TestPBS(unittest.TestCase):

  """Unit tests for PBS integration."""

  def setUp(self):
    """set up a queue to test with."""
    self.qt = _QueueTracker("theq")

  # @mock.patch('subprocess32.check_output')
  @mock.patch('subprocess32.Popen', autospec=True)
  def test_queue_tracker_job_is_running_no_qstat(self, mock_popen):
    """Test that jobs are tracked correctly when no qstat is available."""
    mock_popen.side_effect =\
        build_popen_mock_side_effects({"qsub": ["12547561.hpc-pbs",
                                                "12554695.hpc-pbs"]})

    j1 = PBSJob("j1", "00:20:00", 1, 1, 1000, 1000, 1000, "echo \"test1\"")
    j2 = PBSJob("j2", "00:50:10", 2, 4, 1000, 2000, 5000, "echo \"test2\"")

    self.qt.submit(j1)
    self.qt.submit(j2)
    self.assertEqual(self.qt.is_running(j1.id), True)
    self.assertEqual(self.qt.is_running(j2.id), True)

  @mock.patch('getpass.getuser')
  @mock.patch('time.time')
  @mock.patch('subprocess32.check_output')
  @mock.patch('subprocess32.Popen', autospec=True)
  def test_queue_tracker_job_stat_qs(self, mock_popen, mock_check_output,
                                     time_mock, get_user_mock):
    """Test that jobs are tracked correctly when qstat is available."""
    mock_popen.side_effect =\
        build_popen_mock_side_effects({"qsub": ["12547561.hpc-pbs",
                                                "12554695.hpc-pbs"]})

    qs_head = ("Job ID             Name         User   Time Use S Queue\n" +
               "------------------ ------------ ------ -------- - -----\n")
    qs_j1 = "12547561.hpc-pbs    j1           user1  00:00:01 R theq2\n"
    qs_j2 = "12554695.hpc-pbs    j2           user1  00:00:02 R theq2\n"
    qs_wrong = "12554693.hpc-pbs    j3.pbs       user1  214:08:4 R cmb\n" +\
               "12577676.hpc-pbs    ...str_R2.sh user4  159:33:3 R theq2\n"
    j_subm_count = 0
    time_step = 0
    j1_done = False

    get_user_mock.return_value = "user1"

    def check_output_mock_side_effect(*args, **kwargs):
      opt = {0: qs_head + qs_wrong, 1: qs_head + qs_wrong + qs_j1,
             2: qs_head + qs_wrong + qs_j2
             if j1_done else qs_head + qs_wrong + qs_j1 + qs_j2}
      return opt[j_subm_count]
    mock_check_output.side_effect = check_output_mock_side_effect

    def time_mock_side_effect(*args, **kwargs):
      return {0: 20, 1: 60, 2: 500}[time_step]
    time_mock.side_effect = time_mock_side_effect

    j1 = PBSJob("j1", "00:20:00", 1, 1, 1000, 1000, 1000, "echo \"test1\"")
    j2 = PBSJob("j2", "00:50:10", 2, 4, 1000, 2000, 5000, "echo \"test2\"")

    qt_with_qstat = _QueueTracker("theq2", allow_qstat=True, loud_qstat=False)
    qt_with_qstat.submit(j1)
    j_subm_count += 1
    qt_with_qstat.submit(j2)
    j_subm_count += 1
    # both submitted and running...
    self.assertEqual(qt_with_qstat.is_running(j1.id), True)
    self.assertEqual(qt_with_qstat.is_running(j2.id), True)
    j1_done = True
    time_step += 1
    # j1 finished, but insufficient time passed for checking qstat again.
    self.assertEqual(qt_with_qstat.is_running(j1.id), True)
    self.assertEqual(qt_with_qstat.is_running(j2.id), True)
    # j1 finished, but insufficient time passed for checking qstat again.
    time_step += 1
    self.assertEqual(qt_with_qstat.is_running(j1.id), False)
    self.assertEqual(qt_with_qstat.is_running(j2.id), True)


###############################################################################
#                            EXTERNAL INTERFACE                               #
###############################################################################

# When this module is imported, we instantiate the _PBSTracker class
# to make the following object -- this is the public  itnerface to the module.
pbs_tracker = _PBSTracker()

###############################################################################
#               ENTRY POINT WHEN RUN AS A STAND-ALONE MODULE                  #
###############################################################################

if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
