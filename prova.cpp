#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm> // per std::find
#include <cstdlib>   // per rand() e srand()
#include <ctime>     // per time()
#include <utility>   // per std::pair

// Funzione di crossover tra due array di interi mantenendo l'ordine, escludendo l'ultimo elemento
std::pair<std::vector<int>, std::vector<int>>
orderedCrossover(const std::vector<int> &parent1, const std::vector<int> &parent2) {
    int len = parent1.size();

    // Assicurati che i genitori abbiano la stessa lunghezza e almeno due elementi (per poter escludere l'ultimo)
    if (len != parent2.size() || len < 2) {
        std::cerr << "Parents must have the same length and at least two elements" << std::endl;
        return {};
    }

    // Inizializza la generazione di numeri casuali
    std::srand(std::time(0));

    // Scegli un punto di crossover casuale, escludendo l'ultimo elemento
    int crossover_point = std::rand() % (len - 1);

    // Array per i figli, escluso l'ultimo elemento
    std::vector<int> offspring1(parent1.begin(), parent1.begin() + crossover_point);
    std::vector<int> offspring2(parent2.begin(), parent2.begin() + crossover_point);

    // Crea un insieme per tenere traccia degli elementi già presenti
    std::unordered_set<int> offspring1_set(offspring1.begin(), offspring1.end());
    std::unordered_set<int> offspring2_set(offspring2.begin(), offspring2.end());

    // Aggiungi gli elementi mancanti da parent2 a offspring1 mantenendo l'ordine, escludendo l'ultimo elemento
    for (int i = 0; i < len - 1; ++i) {
        int elem = parent2[i];
        if (offspring1_set.find(elem) == offspring1_set.end()) {
            offspring1.push_back(elem);
            offspring1_set.insert(elem);
        }
    }

    // Aggiungi gli elementi mancanti da parent1 a offspring2 mantenendo l'ordine, escludendo l'ultimo elemento
    for (int i = 0; i < len - 1; ++i) {
        int elem = parent1[i];
        if (offspring2_set.find(elem) == offspring2_set.end()) {
            offspring2.push_back(elem);
            offspring2_set.insert(elem);
        }
    }

    // Aggiungi l'ultimo elemento invariato ai figli
    offspring1.push_back(parent1.back());
    offspring2.push_back(parent2.back());

    return {offspring1, offspring2};
}

int main() {
    // Inizializza due array di genitori
    std::vector<int> parent1 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1};
    std::vector<int> parent2 = {1, 2, 3, 5, 8, 7, 6, 9, 10, 4, 1};

    // Esegui il crossover ordinato
    auto [offspring1, offspring2] = orderedCrossover(parent1, parent2);

    // Stampa i risultati
    std::cout << "Parent 1: ";
    for (int val: parent1) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Parent 2: ";
    for (int val: parent2) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Offspring 1: ";
    for (int val: offspring1) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Offspring 2: ";
    for (int val: offspring2) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}