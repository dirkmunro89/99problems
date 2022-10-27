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
	rm -f *.log
	rm -f *.npz
	rm -f *.png
	rm -f *.tex
	rm -f *.aux
	rm -f *.dat
	rm -f *.eps
	rm -f *.pdf
	rm -f *.pgf
	rm -f *.vtp

.PHONY: docsclean
docsclean:
	rm -rf docs/build/

.PHONY: distclean
distclean:
	rm -rf $(VENV)/
	rm -rf .tox/
	rm -rf isct.egg-info/
	rm -rf cov_html/
