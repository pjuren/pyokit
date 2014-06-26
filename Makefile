




docs :
	$(MAKE) -C sphinx clean
	$(MAKE) -C sphinx html 
	rm -rf docs/doctrees
	mv docs/html/* docs
	mv docs/html/.buildinfo docs
	rmdir docs/html
.PHONY : docs
