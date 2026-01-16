valgrind \
  --tool=memcheck \
  --track-origins=yes \
  --leak-check=full \
  --num-callers=30 \
  ./$1 -s