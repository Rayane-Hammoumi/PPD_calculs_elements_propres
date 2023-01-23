SOURCES = main.c fct_auxiliaires.c
OUTPUT = calcul_elements_propres
LIB = -lm -lgsl -lgslcblas 
OPT = -Wall -g
OBJ = $(SOURCES:.c=.o)

all: $(OUTPUT)

clean: 
	rm $(OUTPUT) *.o ./*.o

$(OUTPUT): $(OBJ)
	@gcc $^  -o $(OUTPUT)  $(LIB)
	rm *.o

%.o: %.c
	@gcc $(OPT) -c $<  -o $@






