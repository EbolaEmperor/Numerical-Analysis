all:
	xelatex solutions

release:
	mv solutions.pdf ..

clean:
	$(shell if [ -e *.bcf ];then rm *.bcf; fi)
	$(shell if [ -e *.log ];then rm *.log; fi)
	$(shell if [ -e *.out ];then rm *.out; fi)
	$(shell if [ -e *.toc ];then rm *.toc; fi)
	$(shell if [ -e *.run.xml ];then rm *.run.xml; fi)
	rm *.aux