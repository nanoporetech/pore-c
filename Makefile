.PHONY: init test clean bumpversion

BUMP_TYPE := micro

test:
	@tox -e py37

clean:
	rm -f MANIFEST
	rm -rf build dist src/*.egg_info .pytest_cache
	find . -name '__pycache__' -exec rm -rf '{}' \;

conda-package:
	conda build conda.recipe

bumpversion:
	@tox -e dephell -- project bump  $(BUMP_TYPE)

package:
	@tox -e dephell -- project build --from=setup.py

