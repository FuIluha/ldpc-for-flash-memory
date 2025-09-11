#include <stdio.h>
#include "progress_bar.h"

void print_progress_bar(size_t count, size_t max)
{
    const char prefix[] = "Progress: [";
    const size_t bar_width = 50;

    size_t filled = (count * bar_width) / max;
    printf("\r%s", prefix);
    for (size_t i = 0; i < bar_width; ++i) {
        putchar(i < filled ? '#' : ' ');
    }
    printf("]%3zu%%", (count * 100) / max);
    if (count == max) {
        printf("\nTask done!\n");
    }
    fflush(stdout);
}