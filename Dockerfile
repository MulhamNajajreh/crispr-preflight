# Lightweight zero-dependency container for Hugging Face Spaces (sdk: docker, port 7860)
FROM python:3.12-slim

WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src ./src

RUN pip install --no-cache-dir . \
    && useradd -m -u 1000 appuser \
    && mkdir -p /home/appuser/.cache/crispr-preflight \
    && chown -R appuser:appuser /home/appuser /app

USER appuser
ENV HOME=/home/appuser
EXPOSE 7860

CMD ["crispr-preflight", "serve", "--host", "0.0.0.0", "--port", "7860"]
