#!/bin/bash

if ! pgrep -x "redis-server" > /dev/null
then
    echo "Redis-server wordt gestart..."
    redis-server --daemonize yes
else
    echo "Redis-server draait al."
fi

echo "Start Celery worker..."
celery -A task worker --loglevel=info &

echo "Start Flask server..."
flask run --host=0.0.0.0 --port=5001 --debug &