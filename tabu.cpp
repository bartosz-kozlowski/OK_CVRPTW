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
#include <utility>

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

bool CheckConditions(const vector<Customer>& CustomersVector, const vector<vector<double>>& distanceMatrix, const vector<Customer>& tempTruckRoute, double& tempSolValue) {
    double AvailableTime = CustomersVector[0].L_end_t;
    double CurrentTime = 0.0;
    int previousPosition = 0;
    for (int i = 1; i < tempTruckRoute.size(); i++)
    {
        int index = tempTruckRoute[i].i;
        if (CurrentTime + distanceMatrix[previousPosition][index] <= CustomersVector[index].L_end_t
            && max(CurrentTime + distanceMatrix[previousPosition][index], (double)CustomersVector[index].e_ready_t) + //Jezeli czas dotarcia/oczekiwania + obslugi + powrotu > mozliwy 
            CustomersVector[index].d_service + distanceMatrix[index][0] <= AvailableTime)
        {
            CurrentTime += max(distanceMatrix[previousPosition][index],
                (double)CustomersVector[index].e_ready_t - CurrentTime) +
                CustomersVector[index].d_service; 
            previousPosition = index;
            continue;
        }
        else {
            return false;
        }
    }

    tempSolValue = CurrentTime + distanceMatrix[previousPosition][0];
    return true;
}

void CalculateTruckValue(const vector<Customer>& CustomersVector, const vector<vector<double>>& distanceMatrix, const vector<Customer>& tempTruckRoute, double& SolValue) {
    double CurrentTime = 0.0;
    int previousPosition = 0;
    for (int i = 1; i < tempTruckRoute.size(); i++)
    {
        int index = tempTruckRoute[i].i;
        
        CurrentTime += max(distanceMatrix[previousPosition][index],
        (double)CustomersVector[index].e_ready_t - CurrentTime) +
        CustomersVector[index].d_service; 
        previousPosition = index;
    }
    SolValue = CurrentTime + distanceMatrix[previousPosition][0];
}

bool Neighbour(const vector<Customer>& CustomersVector, vector<vector<Customer>>& NeighbourSolution, const int& capacity, const vector<vector<double>>& distanceMatrix,
    const vector<vector<Customer>>& SolutionVector, const int& TruckId, const int& CustomerNumInTruckRoute) {

    bool deletedTruck = false;
    int ClientIndex = SolutionVector[TruckId][CustomerNumInTruckRoute].i; // indeks klienta z trasy ciezorwki (CustomerNumInTruckRoute to numer w kolejnoci dla trasy)
    NeighbourSolution = SolutionVector;
    if (NeighbourSolution[TruckId].size() == 2) { // Jezeli size = 2 (magazyn + 1 klient, usuwa ciezarowe)
        NeighbourSolution.erase(NeighbourSolution.begin() + TruckId); //usuwa ciezarowe o id TruckId
        deletedTruck = true;
    }
    else {
        NeighbourSolution[TruckId].erase(NeighbourSolution[TruckId].begin() + CustomerNumInTruckRoute); //usuwa klienta o id CustomerNumInTruckRoute
    }

    double AvailableTime = CustomersVector[0].L_end_t; //czas zamkniecia magazynu
    bool insertSucces = false;

    vector<Customer> bestSolutionVector;
    double bestValue = numeric_limits<double>::max();
    int bestTruckNum;

    for (int i = 0; i < NeighbourSolution.size(); i++) { //ciezarowy

        if (i == TruckId && deletedTruck == false)
        {
            continue;
        }
        int restCapacity = capacity - CustomersVector[ClientIndex].q_demand;
        for (int j = 1; j < NeighbourSolution[i].size(); j++)
        {
            restCapacity -= NeighbourSolution[i][j].q_demand;
        }
        if (restCapacity < 0)
        {
            continue;
        }
        double originalTruckValue;
        CalculateTruckValue(CustomersVector, distanceMatrix, NeighbourSolution[i], originalTruckValue);
        vector<Customer> tempTruckRoute;
        for (int j = 1; j < NeighbourSolution[i].size() + 1; j++) // testujemy po kolei wszystkie miejsca w ciezarowie (od 1 do ostatniego dodatkowego na samym koncu)
        {
            tempTruckRoute.clear();
            tempTruckRoute.push_back(CustomersVector[0]);
            for (int k = 1; k < j; k++) // przepisanie klientow przed miesjscem wsadzenia
            {
                tempTruckRoute.push_back(NeighbourSolution[i][k]);
            }
            tempTruckRoute.push_back(CustomersVector[ClientIndex]); //wsadzenie

            for (int n = j; n < NeighbourSolution[i].size(); n++)//wsadzanie pozostalych
            {
                int index = NeighbourSolution[i][n].i;
                tempTruckRoute.push_back(CustomersVector[index]);
            }
            double tempSolValue = 0.0;
            bool currentInsertSucces;
            currentInsertSucces = CheckConditions(CustomersVector, distanceMatrix, tempTruckRoute, tempSolValue);
            if (currentInsertSucces)
            {
                insertSucces = true;
                if (tempSolValue - originalTruckValue < bestValue)
                {
                    bestValue = tempSolValue - originalTruckValue;
                    bestTruckNum = i;
                    bestSolutionVector = tempTruckRoute;
                }
            }
        }
    }

    if (insertSucces)
    {
        NeighbourSolution[bestTruckNum].clear();
        NeighbourSolution[bestTruckNum] = bestSolutionVector;
        return true;
    }
    return false;
}


void Greedy(const vector<Customer>& CustomersVector, vector<vector<Customer>>& TrucksAndRoutes, const int& capacity, const vector<vector<double>>& distanceMatrix) {
    vector<int>UnvisitedCustomers(CustomersVector.size() - 1);

    for (int i = 1; i < CustomersVector.size(); i++) {
        UnvisitedCustomers[i - 1] = i;
    }

    double AvailableTime = CustomersVector[0].L_end_t;
    int NextCustomerToVisit;
    while (UnvisitedCustomers.size() > 0) {
        int AvailableCapacity = capacity;
        int CurrentPosition = 0;
        vector<Customer> TruckRoute;
        TruckRoute.push_back(CustomersVector[0]);
        double CurrentTime = 0.0;
        bool possibility_to_find = true;
        while (UnvisitedCustomers.size() > 0) {
            double BestValue = numeric_limits<double>::max();
            possibility_to_find = false;

            for (int i = 0; i < UnvisitedCustomers.size(); i++) {
                int CandidateToVisit = UnvisitedCustomers[i];

                double PossibleValue = ((max(distanceMatrix[CurrentPosition][CandidateToVisit],
                    (double)CustomersVector[CandidateToVisit].e_ready_t - CurrentTime) +
                    CustomersVector[CandidateToVisit].d_service)); //Mozliwa wartosc = max ( droga do punktu ALBO droga + czas oczekiwania do rozpoczecia ) + czas obslugi

                if (PossibleValue < BestValue && //Czy jest lepszy od juz najlepszego i czy spelnia wymogi
                    CustomersVector[CandidateToVisit].q_demand <= AvailableCapacity &&
                    CurrentTime + distanceMatrix[CurrentPosition][CandidateToVisit] <= CustomersVector[CandidateToVisit].L_end_t &&
                    max(CurrentTime + distanceMatrix[CurrentPosition][CandidateToVisit], (double)CustomersVector[CandidateToVisit].e_ready_t) +
                    CustomersVector[CandidateToVisit].d_service + distanceMatrix[CandidateToVisit][0] <= AvailableTime) {
                    BestValue = PossibleValue;
                    NextCustomerToVisit = CandidateToVisit;
                    possibility_to_find = true;
                }
            }

            if (possibility_to_find == true) {
                TruckRoute.push_back(CustomersVector[NextCustomerToVisit]);
                UnvisitedCustomers.erase(remove(UnvisitedCustomers.begin(), UnvisitedCustomers.end(), NextCustomerToVisit), UnvisitedCustomers.end());
                //remove - funkcja przesuwa wszystkie elementy o wartosci rownej index na koniec wektora i zwraca index pierwszego z nich (u nas jest tylko jeden taki).
                //erase - funkcja usuwa elementy od miejsca wskazywanego przez remove do konca wektora przez co element o wartosci index zostaje usuniety.
                CurrentTime += max(distanceMatrix[CurrentPosition][NextCustomerToVisit],
                    (double)CustomersVector[NextCustomerToVisit].e_ready_t - CurrentTime) +
                    CustomersVector[NextCustomerToVisit].d_service;
                AvailableCapacity -= CustomersVector[NextCustomerToVisit].q_demand;
                CurrentPosition = NextCustomerToVisit;
            }
            else {
                break;
            }
        }

        TrucksAndRoutes.push_back(TruckRoute); //dodawana jest ciezarowka
    }
}

bool compareSolutions(const vector<vector<Customer>>& sol1, const vector<vector<Customer>>& sol2) {
    if (sol1.size() != sol2.size()) {
        return false;
    }
    for (int i = 0; i < sol1.size(); i++) {
        if (sol1[i].size() != sol2[i].size()) {
            return false;
        }

        for (int j = 0; j < sol1[i].size(); j++) {
            if (!(sol1[i][j].i == sol2[i][j].i)) {
                return false;
            }
        }
    }
    return true;
}
bool isInTabu(const vector<vector<vector<Customer>>>& tabu, const vector<vector<Customer>>& sol) {
    for (int i = 0; i < tabu.size(); i++)
    {
        if (compareSolutions(tabu[i], sol))
        {
            return true;
        }
    }
    return false;
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

            vector<vector<Customer>> BestNeighbourSolVector;
            double Best;
            Greedy(CustomersVector, BestNeighbourSolVector, capacity, distance_matrix); // Wypelnia BestNeighbourSolVector z ktorego startujemy
            Best = SolutionValue(BestNeighbourSolVector); //Sprawdza laczny czas wlasnie znalezionego rozwiazania

            vector<vector<Customer>> NeighbourSolVector; // obliczane rozwiazanie sasiednie
            vector<vector<Customer>> LookingForNeighbours; // rozw z ktorego szukamy siadow
            vector<vector<Customer>> CurrBestNeighbourSolVector = BestNeighbourSolVector; // najlepszy z gorszych w danej iteracji (na poczatku ustawiony na wynik greedy)
            double NBValue; //wartosc rozwiazania sasiada
            bool newBest = true;
            bool isNeihPoss = false;

            vector<vector<vector<Customer>>> tabuList;
            int tabuSize = 30;
            while (true)
            {
                currentTime = high_resolution_clock::now(); //obecny czas
                timeSpan = duration_cast<duration<double>>(currentTime - startTime);// oblicza roznice czasow

                if (timeSpan.count() >= 175.0) {
                    break; // Jezeli czas przekroczony przerywa
                }

                double CurrentIterationBest = numeric_limits<double>::max(); // najlepsza wartosc sasiada w iteracji

                if (newBest == true) { // jezeli jest nowy najlepszy:
                    LookingForNeighbours = BestNeighbourSolVector; // szukamy z najlepszego
                    newBest = false;
                }
                else { // jezeli nie ma nowego najlepszego
                    LookingForNeighbours = CurrBestNeighbourSolVector; // szukamy z najlepszy z gorszych
                }
                
                if (!isInTabu(tabuList, LookingForNeighbours))
                {
                    if (tabuList.size() > tabuSize)
                    {
                        tabuList.erase(tabuList.begin() + 0);
                    }
                    tabuList.push_back(LookingForNeighbours);
                }


                for (int i = 0; i < 60; i++) // parametr - liczba prob poszukiwania sasiadow
                {
                    int TruckId; // id losowanej ciezarowy
                    int CustomerNumInTruckRoute;  // numer w kolejnosci odwiedzanego klienta przez ciezarowke TruckId

                    random_device rd; //losowanie 
                    mt19937 gen(rd());
                    uniform_int_distribution<> distribTruck(0, LookingForNeighbours.size() - 1); // zakres to wszystkie ciezarowki w rozwiazaniu z ktorego szukamy sasiadow
                    TruckId = distribTruck(gen); //losowanie ciezarowy
                    uniform_int_distribution<> distribTruckClient(1, LookingForNeighbours[TruckId].size() - 1); // zakres to wszyscy klienci z ciezarowki (oprocz magazynu)
                    CustomerNumInTruckRoute = distribTruckClient(gen); //losowanie klienta

                    NeighbourSolVector.clear(); // czyszcenie wektora z sasiadem
                    isNeihPoss = Neighbour(CustomersVector, NeighbourSolVector, capacity, distance_matrix, LookingForNeighbours, TruckId, CustomerNumInTruckRoute); // funkcja Neighbour zwraca czy sasiad z wylosowanych danych jest mozliwy oraz go szuka

                    if (isNeihPoss) // jezeli dany sasiad jest mozliwy
                    {
                        NBValue = SolutionValue(NeighbourSolVector); // obliczamy wartosc kosztu dla sasiada

                        if (NBValue < Best) { //jezeli wynik sasiada jest lepszy 
                            Best = NBValue;
                            BestNeighbourSolVector = NeighbourSolVector;
                            newBest = true;
                        }
                        if (NBValue < CurrentIterationBest) {
                           if (!isInTabu(tabuList, NeighbourSolVector))
                           {
                               CurrentIterationBest = NBValue;
                               CurrBestNeighbourSolVector = NeighbourSolVector;
                           }
                        }
                    }
                }
            }

            output << BestNeighbourSolVector.size() << " " << fixed << setprecision(5) << Best << "\n"; //Zapis do pliku wyjsciowego
            for (int i = 0; i < BestNeighbourSolVector.size(); i++)
            {
                for (int j = 1; j < BestNeighbourSolVector[i].size(); j++) //pomijamy magazyn w wyniku
                {
                    output << BestNeighbourSolVector[i][j].i << " ";
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