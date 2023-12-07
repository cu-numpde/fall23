build:
	. VENV/bin/activate && jupyter-book build .
	sed -i 's|"img/|"../img/|g' _build/html/slides/*.html

dev:
	. VENV/bin/activate && jupyter-nbclassic --notebook-dir=. slides/

publish:
	. VENV/bin/activate && ghp-import -n -p -f _build/html

clean:
	rm -r _build
