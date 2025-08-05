# Папка для сборки
BUILD_DIR = ../build

# Исходники (относительно текущего каталога)
SRCS = density_evolution.c progress_bar.c density_evolution_poly.c

# Генерация объектников в build/
OBJS = $(patsubst %.c,$(BUILD_DIR)/%.o,$(SRCS))

# Компилятор
CC = clang

# Флаги компиляции и линковки
CFLAGS = -fPIC -fvisibility=hidden -Wall -O2
LDFLAGS = -shared

# Цель по умолчанию
all: $(BUILD_DIR)/libdensity.so

# Сборка .so
$(BUILD_DIR)/libdensity.so: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

# Создание объектников в build/
$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Убедиться, что папка build существует
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Очистка
clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/libdensity.so