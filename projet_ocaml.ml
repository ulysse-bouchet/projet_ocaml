(* Imports *)
open Format

(* Séquences génétiques - Types de bases *)
type nucleotide = A | C | G | T;;
type brin = nucleotide list;;

(*** Question 1 ***)

(* 
@function contenu_gc 
@desc Pour un brin d'ADN donné, renvoie le taux de nucléotides de type G ou C.
@param b Un brin d'ADN
@returns float Le taux de nucléotides de type G ou C du brin 
*)
let contenu_gc (b: brin): float =
    if List.length b = 0 then 
        failwith "[contenu_gc] erreur : brin vide"
    else
        List.fold_left (fun acc nucl -> acc +. if nucl = G || nucl = C then 1. else 0.) 0. b
        /. 
        float_of_int (List.length b)
    ;;

(* Assertions *)
let () = assert (contenu_gc [A;T;G;T;T;G;A;C] = 0.375);;
let () = assert (contenu_gc [C;T;T;A] = 0.25);;
let () = assert (contenu_gc [A;A;A;T;A] = 0.);;
(* contenu_gc [] provoquerait une erreur car la liste est vide *)
let () = printf "%-30s %s\n" "contenu_gc:" "Assertions effectuées avec succès.";;

(*** Question 2 ***)
(*
@function nucleotide_complementaire
@desc Pour un nucléotide donné, renvoie le nucléotide complémentaire (A <-> T, C <-> G).
@param nuc Un nucléotide 
@returns Le nucléotide complémentaire (A <-> T, C <-> G)
*)
let nucleotide_complementaire (nuc: nucleotide): nucleotide = 
    match nuc with
        | T -> A
        | A -> T
        | G -> C
        | C -> G
    ;;

(* Assertions *)
let () = assert (nucleotide_complementaire T = A);;
let () = assert (nucleotide_complementaire A = T);;
let () = assert (nucleotide_complementaire C = G);;
let () = assert (nucleotide_complementaire G = C);;
let () = printf "%-30s %s\n" "nucleotide_complementaire:" "Assertions effectuées avec succès.";;

(*
@function brin_complementaire
@desc Pour un brin d'ADN donné, renvoie le brin complémentaire (càd. chaque nucléotide est complémentaire).
@param b Un brin d'ADN
@returns Le brin dont chaque nucléotide est complémentaire à celui associé du brin b
*)
let brin_complementaire (b: brin): brin = 
        List.fold_left (fun acc nucl -> acc@[nucleotide_complementaire nucl]) [] b
    ;;

(* Assertions *)
let () = assert (brin_complementaire [T] = [A]);;
let () = assert (brin_complementaire [C;T;T;C] = [G;A;A;G]);;
let () = assert (brin_complementaire [C;T;A;A;T;G;T] = [G;A;T;T;A;C;A]);;
let () = assert (brin_complementaire [] = []);;
let () = printf "%-30s %s\n" "brin_complementaire:" "Assertions effectuées avec succès.";;

(*** Question 3 ***)

(*
@function distance
@desc Pour deux brins donnés, compte le nombre de nucléotides différents.
@param b1 Un brin d'ADN
@param b2 Un autre brin d'ADN
@returns Le nombre de nucléotides différents
@throws Une exception lorsque les brins sont de tailles différentes
*)
let distance (b1: brin) (b2: brin): int = 
    let rec distance_rec (b1: brin) (b2: brin) (cpt: int): int = 
        match (b1, b2) with
        | [], []                -> cpt
        | [], _                 -> failwith "[distance_rec] erreur : brins de taille différente"
        |  _, []                -> failwith "[distance_rec] erreur : brins de taille différente"
        | (h1::t1), (h2::t2)    -> if h1 <> h2 then distance_rec t1 t2 (cpt + 1) else distance_rec t1 t2 cpt
    in
        distance_rec b1 b2 0
    ;;

(* Assertions *)
let () = assert (distance [T] [T] = 0);;
let () = assert (distance [T] [C] = 1);;
let () = assert (distance [G;A;G] [A;G;G] = 2);;
let () = assert (distance [] [] = 0);;
(* distance [A] [A;A] provoquerait une erreur car les listes sont de longueur différentes*)
let () = printf "%-30s %s\n" "distance:" "Assertions effectuées avec succès.";;

