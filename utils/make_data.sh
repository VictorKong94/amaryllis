for FILE in $(ls original); do
    python3 truncate.py original/$FILE 50000
done
