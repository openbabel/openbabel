openbabel:
	cd src
	make -f bc32.mak openbabel

lib:
	cd src
	make -f bc32.mak lib
default: openbabel

tidy:
	cd src
	make -f bc32.mak tidy

clean:
	cd src
	make -f bc32.mak clean

