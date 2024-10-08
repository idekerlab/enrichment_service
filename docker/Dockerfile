FROM continuumio/miniconda3


RUN mkdir /tmp/enrichment_service
COPY dist/*.whl /tmp/enrichment_service/

RUN pip install /tmp/enrichment_service/enrichment_service*whl

RUN rm -rf /tmp/enrichment_service

ENTRYPOINT ["/opt/conda/bin/enrichment_servicecmd.py"]
CMD ["--help"]
