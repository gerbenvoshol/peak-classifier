all:
	gcc -O2 -std=gnu99 libxtend.c biolibc.c peak-classifier.c -o peak-classifier
	gcc -O2 -std=gnu99 libxtend.c biolibc.c filter-overlaps.c -o filter-overlaps

clean:
	$(RM) peak-classifier
	$(RM) filter-overlaps

install:
	cp peak-classifier /usr/local/bin/
	cp filter-overlaps /usr/local/bin/
	cp feature-view.py /usr/local/bin/
	cp extract-genes.* /usr/local/bin/
	cp Man/* /usr/local/man/man1/