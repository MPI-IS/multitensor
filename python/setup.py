import os

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext


try:
    from Cython.Build import cythonize
    USE_CYTHON = True
    print("Found Cython")
except:
    # XXX: When trying to install from the pre-transpiled C sources I
    # get the following warning:
    # >>> import multitensor
    # Traceback (most recent call last):
    #   File "<stdin>", line 1, in <module>
    #   File "/home/me/Code/multi-tensor-factorization/python/venv/lib/python3.8/site-packages/multitensor-1.0.0-py3.8-linux-x86_64.egg/multitensor/__init__.py", line 3, in <module>
    #     from .main import py_multitensor_factorization, PyReport as Report
    # ImportError: dynamic module does not define module export function (PyInit_main)
    USE_CYTHON = False


def get_version():
    """Read the package version from version.txt"""

    version_path = "../version.txt"
    if not os.path.exists(version_path):
        return

    with open(version_path) as fd:
        return fd.read().strip()


class build_ext(_build_ext):
    def finalize_options(self):
        super().finalize_options()
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False

        import numpy
        self.include_dirs.append(numpy.get_include())


EXTRA_COMPILE_ARGS = ["-std=c++17"]

INCLUDE_DIRS = [
    "../applications/include/",
    "../include/"
]
LIBRARY_DIRS = []

# Boost headers
BOOST_INCLUDE_DIR = os.getenv("BOOST_INCLUDE_DIR")
if BOOST_INCLUDE_DIR:
    print("Adding BOOST_INCLUDE_DIR = {} to include dirs".format(BOOST_INCLUDE_DIR))
    INCLUDE_DIRS.append(BOOST_INCLUDE_DIR)

# Boost libraries
BOOST_LIBRARY_DIR = os.getenv("BOOST_LIBRARY_DIR")
if BOOST_LIBRARY_DIR:
    print("Adding BOOST_LIBRARY_DIR = {} to the library dirs".format(BOOST_LIBRARY_DIR))
    LIBRARY_DIRS.append(BOOST_LIBRARY_DIR)


ext = ".pyx" if USE_CYTHON else ".cpp"
extensions = [
    Extension(
        "multitensor.app_utils", [
            "multitensor/app_utils" + ext,
            "../applications/src/app_utils.cpp"
        ],
        extra_compile_args=EXTRA_COMPILE_ARGS,
        include_dirs=INCLUDE_DIRS,
        libraries=["boost_filesystem", "boost_system"],
        library_dirs=LIBRARY_DIRS,
        runtime_library_dirs=LIBRARY_DIRS,
        language="c++",
    ),
    Extension(
        "multitensor.main", [
            "multitensor/main" + ext
        ],
        extra_compile_args=EXTRA_COMPILE_ARGS,
        include_dirs=INCLUDE_DIRS,
        libraries=["boost_graph"],
        library_dirs=LIBRARY_DIRS,
        runtime_library_dirs=LIBRARY_DIRS,
        language="c++",
    ),
]


if USE_CYTHON:
    print("Using Cython to regenerate cpp files.")
    extensions = cythonize(
        extensions,
        compiler_directives={
            "language_level": 3
        }
    )


setup(
    name="multitensor",
    version=get_version(),
    description=(
        "A library for multilayer network tensor factorization that "
        "can be used for community detection, link prediction and measure "
        "layer interdependence."
    ),
    long_description=open("../README.md").read(),
    long_description_content_type="text/markdown",
    url="https://code.is.localnet/project/21/",
    author="Caterina De Bacco, Jean-Claude Passy, Ivan Oreshnikov",
    author_email="caterina.debacco@tuebingen.mpg.de, jean-claude.passy@tuebignen.mpg.de, ivan.oreshnikov@tuebingen.mpg.de",

    python_requires=">=3.5",
    install_requires=[
        "numpy"
    ],
    cmdclass={  # this is necessary to find numpy headers directory.
        "build_ext": build_ext
    },
    ext_modules=extensions,

    packages=["multitensor"],
)
