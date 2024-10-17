#!/bin/bash

echo "--------------------Uninistall any old versions...-----------------"
pip3 uninstall -y pycheron

echo "--------------------Building Pycheron----------------------"
poetry build

echo "--------------------Installing Pycheron--------------------"
poetry run pip3 install --upgrade pip
poetry run pip3 install numpy
poetry run pip3 install llvmlite
poetry run pip3 install ./dist/pycheron-3.0.0-py3-none-any.whl
