## Documentation is available at

https://mstekiel.github.io/crysfipy/

It is formatted according to the numpydoc standards and Google Style for Python

https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard

https://djangocas.dev/docs/4.0/python-docstring-google-style-example.html

### Compiling and setting up the documentation.

For changes check the `index.rst` file and the contents od the `_build` folder.

For windows run in terminal
>./make.bat html

Alternatively, in the environment where sphinx is installed run
> sphinx-build -M html . build

### Development guide

Some useful links:
 - lot of definitions and resources (https://lpn-doc-sphinx-primer.readthedocs.io/en/stable/concepts.html)
 - referencing with sphinx (https://www.sphinx-doc.org/en/master/usage/referencing.html#role-ref)


### Other
Some of the build directories have been rearranged in order to be able to host the documentation on GitHub Pages.
The main thing is the redirecting index.html file in the main folder, that redirects to the real index file.

This idea is from the blog post:
https://python.plainenglish.io/how-to-host-your-sphinx-documentation-on-github-550254f325ae

Based on the thread:
https://github.com/sphinx-doc/sphinx/issues/3382

From which I also incorporated the idea of keeping the build in the build folder, and adding the empty .nojekyll file to prevent GitHub from running it.