




docs :
	$(MAKE) -C sphinx clean
	$(MAKE) -C sphinx html 
	rm -rf Docs/doctrees
	mv Docs/html/* docs
	mv Docs/html/.buildinfo docs
	rmdir Docs/html
.PHONY : docs
