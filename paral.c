#include <pthread.h>

void* func(void *arg)
{
    while (1);
}

int main()
{
    #define NUM_THREADS 8 //use the number of cores (if known)
    pthread_t threads[NUM_THREADS];

    for (int i=0; i < NUM_THREADS; ++i)
        pthread_create(&threads[i], NULL, func, NULL);

    for (int i=0; i < NUM_THREADS; ++i)
        pthread_join(threads[i], NULL);

    return 0;
}