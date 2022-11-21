export VENV := ./venv

.PHONY: install
install: python

.PHONY: python
python:  venv
	. $(VENV)/bin/activate && pip install -e .[dev]

venv:
	test -d $(VENV) || python3 -m venv $(VENV)
	. $(VENV)/bin/activate && pip install --upgrade pip setuptools wheel

.PHONY: test
test:
	tox

.PHONY: docs
docs:
	sphinx-build -b html docs/source/ docs/build/

.PHONY: clean
clean:
	-mv  *.log ./trsh/ 
	-mv  *.npz ./trsh/ 
	-mv  *.png ./trsh/
	-mv  *.tex ./trsh/
	-mv  *.aux ./trsh/
	-mv  *.dat ./trsh/
	-mv  *.eps ./trsh/
	-mv  *.pdf ./trsh/
	-mv  *.pgf ./trsh/
	-mv  *.vtp ./trsh/

.PHONY: docsclean
docsclean:
	rm -rf docs/build/

.PHONY: distclean
distclean:
	rm -rf $(VENV)/
	rm -rf .tox/
	rm -rf isct.egg-info/
	rm -rf cov_html/
