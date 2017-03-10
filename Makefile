
pull:
	git reset --hard origin/master
	git pull

package:
	python3 setup.py sdist

install:
	sudo pip3 install phievo --no-index --find-links file:///home/adrien/Desktop/test-phievo/phievo/dist/phievo-1.0.tar.gz

uninstall:
	sudo pip3 uninstall  phievo
