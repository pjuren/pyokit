




docs :
	$(MAKE) -C sphinx clean
	$(MAKE) -C sphinx html 
	rm -rf Docs/doctrees
	mv Docs/html/* Docs
	mv Docs/html/.buildinfo Docs
	rmdir Docs/html
.PHONY : docs
