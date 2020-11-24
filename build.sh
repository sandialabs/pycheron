#!/bin/bash

version="$1"

if [[ $# -eq 0 ]] ; then
    echo "Error: Must supply a version number"
    exit 1
fi

echo "--------------------Uninistall any old versions...-----------------"
pip3 install --upgrade pip
pip3 uninstall -y pycheron

echo "--------------------Building Pycheron $version----------------------"
poetry build

echo "--------------------Install dephell to update setup.py...----------------------"
pip3 install dephell==0.8.3
dephell deps convert --from=pyproject.toml --to=setup.py

echo "--------------------Installing Pycheron $version--------------------"
poetry run pip3 install numpy==1.19.1
poetry run pip3 install llvmlite==0.34.0
poetry run pip3 install ./dist/pycheron-$version-py3-none-any.whl
