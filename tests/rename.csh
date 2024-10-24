#!/bin/bash

for file in carma_hatp12b*.F90
do
    mv -i "${file}" "${file/hatp12b/bdlt}"
done
