.PHONY: all lint test test-cov install dev clean distclean html serve

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	py.test --cov=q2_fmt

install: all
	$(PYTHON) setup.py install

dev: all
	pip install -e .


distclean: ;

html:
	cd book/q2_fmt_book && q2doc autodoc --plugin fmt --output reference .
	cd book/q2_fmt_book && jupyter book build --html
	cp -r book/q2_fmt_book/data/ book/q2_fmt_book/_build/html/data/

serve:
	npx serve book/q2_fmt_book/_build/html/ -p 4000

clean:
	rm -rf book/q2_fmt_book/_build/html/
