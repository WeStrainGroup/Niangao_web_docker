FROM continuumio/miniconda3

WORKDIR /app

COPY conda-env.yaml /tmp/env.yaml
COPY download_databases.sh /tmp/
RUN conda env create -f /tmp/env.yaml --verbose

SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

RUN /bin/bash /tmp/download_databases.sh

ENV TAXONKIT_DB=/opt/conda/envs/myenv/.taxonkit

COPY ./app/ /app/
RUN chmod +x /app/cap3

CMD ["conda", "run", "-n", "myenv", "Rscript", "/app/app.R"]