# 🚚 OK_CVRPTW — Tabu Search dla problemu marszrutyzacji z oknami czasowymi

Projekt realizuje rozwiązanie **problemu marszrutyzacji pojazdów z ograniczoną pojemnością i oknami czasowymi** (*Capacitated Vehicle Routing Problem with Time Windows – CVRPTW*) z wykorzystaniem metaheurystyki **Tabu Search**.

Został wykonany w ramach projektu akademickiego z optymalizacji kombinatorycznej — nasze podejście pobiło wcześniejsze rekordy uczelniane dla wybranych instancji benchmarkowych.

---

## 📌 Opis problemu

Dane wejściowe:
- Flota pojazdów o **ograniczonej pojemności**
- Zbiór **klientów**, z których każdy posiada:
  - określone **zapotrzebowanie**
  - **przedział czasowy**, w którym należy go obsłużyć
- **Magazyn główny** jako punkt startowy i końcowy tras

🎯 Celem jest wyznaczenie **optymalnych tras dla pojazdów**, które:
- Minimalizują **łączny dystans**
- Spełniają ograniczenia **czasowe** i **pojemnościowe**

---

## 🧠 Zastosowany algorytm

Zaimplementowano algorytm **Tabu Search**, oparty na:
- Heurystyce początkowej (*greedy*)
- Generowaniu sąsiedztwa z uwzględnieniem tabu listy
- Lokalnych poprawkach trasy
- Mechanizmach unikania lokalnych minimów
