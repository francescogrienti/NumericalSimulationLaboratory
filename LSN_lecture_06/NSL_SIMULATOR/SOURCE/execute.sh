#!/bin/bash

max_temperature=2.0

temporary_temperature=0.3

function execute_program() {
    make
    if [ $? -eq 0 ]; then
        echo "Compilazione completata con successo. Eseguo il programma..."

        ./simulator.exe

        if [ $? -eq 0 ]; then
            echo "Il programma è stato eseguito correttamente."
        else
            echo "Si è verificato un errore durante l'esecuzione del programma."
        fi
    else
        echo "Si è verificato un errore durante la compilazione del programma."
    fi
}

while (( $(echo "$temporary_temperature <= $max_temperature" | bc -l) )); do
    execute_program
    temporary_temperature=$(echo "$temporary_temperature + 0.10" | bc)
done
