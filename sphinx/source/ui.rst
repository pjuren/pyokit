====================================
Building user interfaces with Pyokit
====================================
Pyokit contains an extension to the built-in python command-line parser. The
pyokit.interface.cli.CLI class is used to instantiate user interfaces, and once
you have one, you can parse a command line. Options for the user interface are
described by objects created from the pyokit.interface.cli.Option class

Here's an example of creating a simple program with a command-line interface.

::

  import os, sys
  from pyokit.interface.cli import CLI, Option

  def getUI():
    programName = os.path.basename(sys.argv[0])
    longDescription  =  "This script foos the biz baz and outputs a bar.\n" +\
                        "Usage: ./" + programName + " <biz.txt> <baz.txt>"
    shortDescription =  longDescription
    ui = CLI(programName, shortDescription, longDescription)
    ui.minArgs = 2
    ui.maxArgs = 2
    ui.addOption(Option(short = "v", long = "verbose",
                        description = "output status messages to stderr ",
                        required = False))
    ui.addOption(Option(short = "n", long = "num", argName = "fooNum",
                        description = "number of times to foo the biz baz",
                        required = False, type = int))
    ui.addOption(Option(short = "o", long = "output_fn", argName = "filename",
                        description = "output bar to this file, else stdout",
                        required = False, type = str))
    ui.addOption(Option(short="h", long="help",
                        description="show this help message ", special=True))
    ui.parseCommandLine(sys.argv[1:])
    return ui

  def _main():
    ui = getUI()
    if ui.optionIsSet("help") :
      ui.usage()
      sys.exit()
    verbose = (ui.optionIsSet("verbose") == True)
    outfh = open(ui.getValue("output_fn")) \
            if ui.optionIsSet("output_fn") \
            else sys.stdout
    num = ui.getValue("num") if ui.optionIsSet("num") else 1
    biz = [l.strip() for l in open(ui.getArgument(0)) if l.strip() != ""]
    baz = [l.strip() for l in open(ui.getArgument(1)) if l.strip() != ""]
    for i in range(0, num) :
      if (verbose) : sys.stderr.write("fooing the biz baz\n")
      outfh.write("bar " + str(i+1) + ": " +\
                  ",".join(biz) + " -- " + ",".join(baz) + "\n")

  if __name__ == "__main__":
    _main()

.. autoclass:: pyokit.interface.cli.CLI
.. autoclass:: pyokit.interface.cli.Option
