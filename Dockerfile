FROM python:3.7
WORKDIR /app
COPY requirements.txt /app
RUN pip install autopep8 flake8
RUN pip install -r requirements.txt
