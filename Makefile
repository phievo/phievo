PIP=pip3


pull:
	git reset --hard origin/master
	git pull

push:
	git push origin master

package:
	python3 setup.py sdist

install_dependencies:
	sudo $(PIP) install -r requirements.txt

install:
	sudo $(PIP) install phievo --no-index --find-links file:///home/adrien/Desktop/test-phievo/phievo/dist/phievo-1.0.tar.gz

uninstall:
	sudo $(PIP) uninstall  phievo
