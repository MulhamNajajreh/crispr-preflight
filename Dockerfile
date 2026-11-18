# Lightweight zero-dependency container (Render, Cloud Run, Railway, or self-hosted Docker)
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
ENV PORT=7860
EXPOSE 7860

CMD ["sh", "-c", "crispr-preflight serve --host 0.0.0.0 --port ${PORT:-7860}"]
