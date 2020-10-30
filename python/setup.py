"""
Configuration to create a wheel.

:note: got help from
https://stackoverflow.com/questions/35112511/pip-setup-py-bdist-wheel-no-longer-builds-forced-non-pure-wheels/36886459#36886459
"""

from setuptools import setup, Distribution


class BinaryDistribution(Distribution):
    """Distribution which always forces a binary package with platform name."""

    def has_ext_modules(self, *args, **kwargs):
        return True


setup(
    name="multitensor",
    version="{{version}}",
    description=(
        "A library for multilayer network tensor factorization that "
        "can be used for community detection, link prediction and measure "
        "layer interdependence."
    ),
    url="https://code.is.localnet/project/21/",
    author="Caterina De Bacco, Jean-Claude Passy, Ivan Oreshnikov",
    author_email=(
        "caterina.debacco@tuebingen.mpg.de, "
        "jean-claude.passy@tuebignen.mpg.de, "
        "ivan.oreshnikov@tuebingen.mpg.de"
    ),
    python_requires=">=3.6",
    install_requires=["numpy"],
    packages=['.'],
    package_data={".": [
        "{{libname}}"
    ]},
    distclass=BinaryDistribution
)
