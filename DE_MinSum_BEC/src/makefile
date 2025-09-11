# Папка для сборки
BUILD_DIR = ../build

# Исходники
SRCS = density_evolution.c progress_bar.c density_evolution_poly.c differential_evolution_poly_find.c differential_evolution_prefect_border_find.c

# Объектники
OBJS = $(patsubst %.c,$(BUILD_DIR)/%.o,$(SRCS))

# Компилятор
CC = clang

# Пути к OpenMP (brew install libomp)
OMP_PREFIX = /opt/homebrew/opt/libomp
OMP_INC = $(OMP_PREFIX)/include
OMP_LIB = $(OMP_PREFIX)/lib

# Флаги
CFLAGS = -fPIC -fvisibility=hidden -Wall -O2 -Xpreprocessor -fopenmp -I$(OMP_INC)
LDFLAGS = -shared -L$(OMP_LIB) -lomp

# Цель
all: $(BUILD_DIR)/libdensity.so

$(BUILD_DIR)/libdensity.so: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/libdensity.so