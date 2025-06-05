// #include <pthread.h>
// #include <stdio.h>
// #include <stdatomic.h>

// typedef struct {
//     long *value;
//     long thread_id;
//     pthread_mutex_t *lock;
// } ThreadData;

// atomic_long counter = 0;

// void *computation(void *add);

// void *increment(void *arg);

// int main() {
//     pthread_t thread1;
//     pthread_t thread2;

//     long value = 1;

//     pthread_mutex_t lock;
//     pthread_mutex_init(&lock, NULL);

//     ThreadData data1 = {&value, 1, &lock};
//     ThreadData data2 = {&value, 2, &lock};

//     pthread_create(&thread1, NULL, increment, (void *)&data1);
//     pthread_create(&thread2, NULL, increment, (void *)&data2);

//     pthread_join(thread1, NULL);
//     pthread_join(thread2, NULL);

//     printf("Counter: %d\n", counter);

//     pthread_mutex_destroy(&lock);
//     return 0;
// }

// void *computation(void *arg) {
//     ThreadData *data = (ThreadData *)arg;

//     for (long i = 0; i < 1000; i++) {
//         pthread_mutex_lock(data->lock);
//         (*data->value)++;
//         printf("Thread %ld: %ld\n", data->thread_id, (*data->value));
//         pthread_mutex_unlock(data->lock);
//     }
//     return NULL;
// }

// void *increment(void *arg) {
//     for (int i = 0; i < 10000000; i++) {
//         atomic_fetch_add(&counter, 1);
//     }
//     return NULL;
// }

// #include <pthread.h>
// #include <stdio.h>

// int buffer = 0;
// int is_full = 0;

// pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
// pthread_cond_t cond_produce = PTHREAD_COND_INITIALIZER;
// pthread_cond_t cond_consume = PTHREAD_COND_INITIALIZER;

// void *producer(void *arg) {
//     for (int i = 0; i < 1000000000; i++) {
//         pthread_mutex_lock(&mutex);

//         while (is_full)
//             pthread_cond_wait(&cond_produce, &mutex);

//         buffer = i + 1;
//         is_full = 1;
//         printf("Produced: %d\n", buffer);

//         pthread_cond_signal(&cond_consume);
//         pthread_mutex_unlock(&mutex);
//     }
//     return NULL;
// }

// void *consumer(void *arg) {
//     for (int i = 0; i < 1000000000; i++) {
//         pthread_mutex_lock(&mutex);

//         while (!is_full)
//             pthread_cond_wait(&cond_consume, &mutex);

//         printf("Consumed: %d\n", buffer);
//         is_full = 0;

//         pthread_cond_signal(&cond_produce);
//         pthread_mutex_unlock(&mutex);
//     }
//     return NULL;
// }

// int main() {
//     pthread_t t1, t2;
//     pthread_create(&t1, NULL, producer, NULL);
//     pthread_create(&t2, NULL, consumer, NULL);

//     pthread_join(t1, NULL);
//     pthread_join(t2, NULL);
//     return 0;
// }

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/sysinfo.h>

void *burn_cpu(void *arg) {
    while (1) {
        // Waste CPU cycles
        for (volatile int i = 0; i < 1000000; ++i);
    }
    return NULL;
}

int main() {
    int num_cores = get_nprocs();  // or use sysconf(_SC_NPROCESSORS_ONLN)
    printf("Starting %d threads\n", num_cores);

    pthread_t *threads = malloc(num_cores * sizeof(pthread_t));

    for (int i = 0; i < num_cores; i++) {
        pthread_create(&threads[i], NULL, burn_cpu, NULL);
    }

    for (int i = 0; i < num_cores; i++) {
        pthread_join(threads[i], NULL);
    }

    free(threads);
    return 0;
}
