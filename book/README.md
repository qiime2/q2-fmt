# q2-fmt book

User documentation for the q2-fmt project.

## Usage

### Building the book

If you'd like to develop and/or build the q2-fmt book book, you should:

1. Clone this repository
2. Run `pip install -r book/requirements.txt` (it is recommended you do this within a virtual environment)
3. (Optional) Edit the books source files located in the `book/q2_fmt_book/` directory
4. Run `jupyter-book clean book/q2_fmt_book/` to remove any existing builds
5. Run `jupyter-book build book/q2_fmt_book/`

A fully-rendered HTML version of the book will be built in `book/q2_fmt_book/_build/html/`.

## Credits

This project is created using the excellent open source [Jupyter Book project](https://jupyterbook.org/) and the [executablebooks/cookiecutter-jupyter-book template](https://github.com/executablebooks/cookiecutter-jupyter-book).
