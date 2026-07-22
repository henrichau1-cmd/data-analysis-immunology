.PHONY: data process figures test all clean

data:
	python -m src.data.fetch

process:
	python -m src.data.clean

figures:
	python -m src.viz.plots

test:
	pytest -v

all: data process figures

clean:
	rm -rf data/interim/* data/processed/*.csv figures/*.png
