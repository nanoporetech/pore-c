.PHONY: init test clean

init:
	conda env create

reformat:
	black src --target-version py37 --line-length 120
	isort --recursive src

check:
	flake8 src

test:
	python -m pytest -s tests --basetemp=test_output

clean:
	rm -f MANIFEST
	rm -rf build dist src/*.egg_info .pytest_cache
	find . -name '__pycache__' -exec rm -rf '{}' \;

conda-package:
	conda build conda.recipe
