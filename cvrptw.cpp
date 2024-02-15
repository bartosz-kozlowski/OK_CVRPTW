#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <iomanip>
#include <random>
#include <climits>

using namespace std;
using namespace std::chrono;

class Customer {
public:
    int i; //numer wierzcholka (klienta)
    int x; //wspolrzedna x
    int y; //wspolrzedna y
    int q_demand; //zapotrzebowanie na towar
    int e_ready_t; //poczatek okna
    int L_end_t; //koniec okna czasowego
    int d_service; //czas rozladunku
};

void InputData(const string& nazwa_pliku, vector<Customer>& CustomersVector, int& capacity) {
    ifstream plik(nazwa_pliku); //otwarcie pliku
    string smietnik;

    if (!plik.is_open()) { //czy plik otwarty
        cout << "Nie mozna otworzyc pliku z danymi wejsciowymi " << nazwa_pliku << endl;
        exit(0);
        //return;
    }
    for (int i = 0; i < 4; i++) {
        getline(plik, smietnik); //skip bo niepotrzebne
    }

    plik >> smietnik >> capacity; //ograniczenie ciezarowek do kosza, pojemnosc pobieramy

    for (int i = 0; i < 4; i++) {
        getline(plik, smietnik); //skip bo niepotrzebne
    }

    Customer customer;
    //tutaj dodajemy dane
    while (plik >> customer.i >> customer.x >> customer.y >> customer.q_demand >> customer.e_ready_t >> customer.L_end_t >> customer.d_service) {
        CustomersVector.push_back(customer);
    }

    plik.close();
}

double distanceAB(const Customer& A, const Customer& B) { //Liczy dystans miedzy A i B
    double distance = sqrt((((B.x - A.x) * (B.x - A.x)) + ((B.y - A.y) * (B.y - A.y))));
    return distance;
}

void calculateDistanceMatrix(const vector<Customer>& CustomersVector, vector<vector<double>>& distance_mat) { //Tworzy macierz odleglosci z wektorow
    int N = CustomersVector.size();
    distance_mat.assign(N, vector<double>(N)); //zainicjowanie distance_matrix N-toma vektorami o N elemntach typu double

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                distance_mat[i][j] = 0.0;
            }
            else {
                distance_mat[i][j] = distanceAB(CustomersVector[i], CustomersVector[j]);
            }
        }
    }
}

bool IsSolutionPossible(const vector<Customer>& CustomersVector, const int& capacity, const vector<vector<double>>& distanceMatrix) { //Sprawdza czy w ogole mozna stworzyc rozwiazanie
    double RouteToAndBack;
    for (int i = 1; i < CustomersVector.size(); i++) {
        if (CustomersVector[i].q_demand > capacity) {
            return false;
        }
        RouteToAndBack = max(distanceMatrix[0][i], (double)CustomersVector[i].e_ready_t) + CustomersVector[i].d_service + distanceMatrix[i][0];
        if (distanceMatrix[0][i] > CustomersVector[i].L_end_t || RouteToAndBack > CustomersVector[0].L_end_t) { //Czy wszedzie da sie pojechac przed zakonczeniem okna, rozladowac i wrocic do magazynu
            return false;
        }
    }
    return true;
}

double ServiceStartTime(vector<Customer>& Route, int CustomerNum) {
    if (CustomerNum == 0) {
        return 0; // Czas "obslugi" magazynu to poczatek
    }
    double value = max(ServiceStartTime(Route, CustomerNum - 1) + Route[CustomerNum - 1].d_service +
        distanceAB(Route[CustomerNum - 1], Route[CustomerNum]), double(Route[CustomerNum].e_ready_t)); //wzor rekurencyjny z polecenia
        //b_i=max{b_{i-1}+d_{i-1}+c_{i-1,i},e_i}, b_0=0.
    return value;
}

double SolutionValue(vector<vector<Customer>>& TrucksAndRoutes) { //Oblicza wartosc danego rozwiazania
    double totalValue = 0;
    for (int i = 0; i < TrucksAndRoutes.size(); i++) { // sumujemy czas wszystkich ciezarowek
        int LastCustomerFromRoute = TrucksAndRoutes[i].size() - 1;
        double routeValue = ServiceStartTime(TrucksAndRoutes[i], LastCustomerFromRoute) +
            TrucksAndRoutes[i][LastCustomerFromRoute].d_service +
            distanceAB(TrucksAndRoutes[i][LastCustomerFromRoute], TrucksAndRoutes[i][0]); // 
        //routevalue = Moment rozpoczecia ladunku + obsluga ladunku + powrot do magazynu dla ostatniego klienta. 
        //Moment rozpoczenia rozladunku wyznaczany rekurencyjnie za pomoca ServiceStart
        totalValue += routeValue;// do sumy ostatecznej dodajemy czas dla ciezarowki
    }
    return totalValue;
}

void GRASP(const vector<Customer>& CustomersVector, vector<vector<Customer>>& TrucksAndRoutes, const int& capacity, const vector<vector<double>>& distanceMatrix) {
    int parametr = (CustomersVector.size() / 10) + 1; //parametr
    vector<int> UnvisitedCustomers(CustomersVector.size() - 1); //Vector o rozmiarze = ilosc customerow -1 ( bo 1 to magazyn)
    for (int i = 1; i < CustomersVector.size(); i++) {
        UnvisitedCustomers[i - 1] = i; //Wypelnianie UnvisitedCustomers numerami customerow
    }
    vector<int> PickedCustomers; //Vector wylosowanych customerow
    vector<int> AvailableCustomers; //Vector mozliwych do odwiedzenia z danego punktu
    double AvailableTime = CustomersVector[0].L_end_t;  // Czas jaki trzeba dopilnowac zeby wrocic do magazynu
    while (UnvisitedCustomers.size() > 0) { //Jezeli sa jeszcze nieodwiedzeni customersi to szukaj, kazde uruchomienie petli to nowa ciezarowka

        int AvailableCapacity = capacity; // Dostepna pojemnosc, nowa ciezarowka
        int CurrentPosition = 0; // Obecna pozycja

        vector<Customer> TruckRoute;// droga obecnej ciezarowki
        TruckRoute.push_back(CustomersVector[0]);
        double CurrentTime = 0.0; //Obecny czas
        double BestValue;
        double PossibleValue = 0; //Mozliwy koszt pelnej obslugi sprawdzanego klienta
        int NextCustomerToVisit = 0;
        bool possibility_to_find = true; // domyslne (bylo sprawdzane IsSolutionPossible, jestesmy w magazynie)

        while (possibility_to_find && UnvisitedCustomers.size() > 0) {
            possibility_to_find = false; //ponizej w petli sprawdzimy czy mozna sie gdzies dalej dostac
            for (int i = 0; i < UnvisitedCustomers.size(); i++)
            {
                if (CustomersVector[UnvisitedCustomers[i]].q_demand > AvailableCapacity ||
                    CurrentTime + distanceMatrix[CurrentPosition][UnvisitedCustomers[i]] > CustomersVector[UnvisitedCustomers[i]].L_end_t || //Jezeli zapotrz. > dost. pojemn. || czas dojazdu > koniec okna
                    max(CurrentTime + distanceMatrix[CurrentPosition][UnvisitedCustomers[i]], (double)CustomersVector[UnvisitedCustomers[i]].e_ready_t) +
                    CustomersVector[UnvisitedCustomers[i]].d_service + distanceMatrix[UnvisitedCustomers[i]][0] > AvailableTime) {
                    //Jezeli czas dotarcia/oczekiwania + obslugi + powrotu do mag > mozliwy
                    continue;
                }
                else {
                    AvailableCustomers.push_back(UnvisitedCustomers[i]); //Jezeli dany klient spelnia warunki dodajemy go do wektora mozliwych
                    possibility_to_find = true; //mozliwe jest znalezienie nowego klienta
                }
            }

            if (possibility_to_find) {

                if (AvailableCustomers.size() <= parametr) { //Jezeli liczba customersow nieodwiedzonych jest mniejsza lub rowna parametrowi to bierz wszystkie
                    for (int i = 0; i < AvailableCustomers.size(); i++) {
                        PickedCustomers.push_back(AvailableCustomers[i]);
                    }
                }
                else {
                    random_device rd; //losowanie z mozliwych do odwiedzenia //obiekt do generowania pseudolosowych liczb
                    mt19937 gen(rd()); //generator
                    while (PickedCustomers.size() < parametr) {
                        uniform_int_distribution<> distrib(0, AvailableCustomers.size() - 1);
                        int randomCustomer = distrib(gen);

                        if (find(PickedCustomers.begin(), PickedCustomers.end(), AvailableCustomers[randomCustomer]) != PickedCustomers.end()) {
                            continue; //Jezeli wylosowany jest w puli to losuj nowego
                        }
                        else {
                            PickedCustomers.push_back(AvailableCustomers[randomCustomer]); //Jezeli wylosowany nie jest w puli to go dodaj
                        }
                    }
                }
                BestValue = numeric_limits<double>::max(); //Przypisana maksymalna wartosc dla double
                int index; //zmienna pomocnicza
                for (int i = 0; i < PickedCustomers.size(); i++) { //wyznaczenie najlepszego do odwiedzenia
                    index = PickedCustomers[i];
                    PossibleValue = ((max(distanceMatrix[CurrentPosition][index],
                        (double)CustomersVector[index].e_ready_t - CurrentTime) +
                        CustomersVector[index].d_service)); //Mozliwa wartosc = max ( droga do punktu ALBO droga + czas oczekiwania do rozpoczecia ) + czas obslugi
                    if (PossibleValue < BestValue) { //jezeli mozliwy jest lepszy to go ustaw jako najlepszy
                        BestValue = PossibleValue;
                        NextCustomerToVisit = PickedCustomers[i]; // ustaw go jako najlepszy do odwiedzenia/ nastepny
                    }
                }

                index = NextCustomerToVisit;
                TruckRoute.push_back(CustomersVector[index]);
                UnvisitedCustomers.erase(remove(UnvisitedCustomers.begin(), UnvisitedCustomers.end(), index), UnvisitedCustomers.end());
                //remove - funkcja przesuwa wszystkie elementy o wartosci rownej index na koniec wektora i zwraca index pierwszego z nich (u nas jest tylko jeden taki).
                //erase - funkcja usuwa elementy od miejsca wskazywanego przez remove do konca wektora przez co element o wartosci index zostaje usuniety.
                CurrentTime += max(distanceMatrix[CurrentPosition][index],
                    (double)CustomersVector[index].e_ready_t - CurrentTime) +
                    CustomersVector[index].d_service; //obecny czas = po obsluzeniu customera
                AvailableCapacity -= CustomersVector[index].q_demand; // zmniejsz ladunek
                CurrentPosition = index; //Obecna pozycja ustawiona na obsluzonego customera
                PickedCustomers.clear(); //Czyscimy dla nowego klienta
                AvailableCustomers.clear(); //Czyscimy dla nowego klienta
            }
            else {
                break; //Jezeli dana ciezarowka juz nigdzie nie moze pojechac (np. skonczyl sie ladunek lub nie zdazy wrocic do magazynu) to ja zawroc i wyslij nowa
            }
        }
        TrucksAndRoutes.push_back(TruckRoute);
    }
}

int main(int argc, char** argv) {
    //Pomiar czasu od samego poczatku
    high_resolution_clock::time_point startTime = high_resolution_clock::now(); // zapisz czas startu (punkt czasowy o wysokiej dokladnosci/rozdzielczosci)
    high_resolution_clock::time_point currentTime; 
    duration<double> timeSpan; // pozwala na obliczanie roznicy czasow i zapisywaniu liczby sekund w typie double

    vector<Customer> CustomersVector; //wektor przechowujacy klientow
    vector<vector<double>> distance_matrix;
    int capacity;

    InputData(argv[1], CustomersVector, capacity);

    calculateDistanceMatrix(CustomersVector, distance_matrix);//oblicz i wypelnij macierz odleglosci

    bool isPossible = IsSolutionPossible(CustomersVector, capacity, distance_matrix); //Czy rozwiazanie jest mozliwe

    // Zapisz wynik do pliku
    ofstream output("solution_file.txt");

    if (output.is_open()) {
        if (isPossible) { //Jezeli trasa jest dopuszczalna
            vector<vector<Customer>> SolutionVector; //Wektor z roziazaniem
            double TempValue;
            double bestSolutionValue = numeric_limits<double>::max(); // Deklaracja bestSolutionValue, przypisana najwieksza wartosc jaka moze przyjac double

            while (true) {
                currentTime = high_resolution_clock::now();
                timeSpan = duration_cast<duration<double>>(currentTime - startTime);// oblicza roznice czasow

                if (timeSpan.count() >= 170.0) {
                    break; // Jezeli czas przekroczony przerywa
                }

                vector<vector<Customer>> TempSolutionVector;
                GRASP(CustomersVector, TempSolutionVector, capacity, distance_matrix); // Wypelnia TempSolutionVector
                TempValue = SolutionValue(TempSolutionVector); //Sprawdza laczny czas wlasnie znalezionego rozwiazania

                if (TempValue < bestSolutionValue) { //Jezeli znalezione rozwiazanie jest lepsze ustawia jako najlepsze
                    bestSolutionValue = TempValue;
                    SolutionVector = TempSolutionVector;
                }
            }

            output << SolutionVector.size() << " " << fixed << setprecision(5) << bestSolutionValue << "\n"; //Zapis do pliku wyjsciowego
            for (int i = 0; i < SolutionVector.size(); i++)
            {
                for (int j = 1; j < SolutionVector[i].size(); j++) //pomijamy magazyn w wyniku
                {
                    output << SolutionVector[i][j].i << " ";
                }
                output << endl;
            }
        }
        else {
            output << "-1\n";
        }
        output.close();
    }
    else {
        cout << "Blad podczas zapisu wyniku do pliku." << endl;
    }
    return 0;
}