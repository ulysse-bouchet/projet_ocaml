(* Imports *)
open Format

(* Séquences génétiques - Types de bases *)
type nucleotide = A | C | G | T;;
type brin = nucleotide list;;

(*** Question 1 ***)

(* 
@function contenu_gc
@param b Un brin d'ADN
@returns float Le taux de nucléotides de type G ou C du brin 
*)
let contenu_gc (b: brin): float =
    if List.length b = 0 then 
        failwith "erreur : brin vide"
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
let () = printf "%-30s %s\n" "contenu_gc:" "Assertions effectuées avec succès."

(*** Question 2 ***)
(*
@function nucleotide_complementaire
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
let () = printf "%-30s %s\n" "nucleotide_complementaire:" "Assertions effectuées avec succès."

(*
@function brin_complementaire
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
let () = printf "%-30s %s\n" "brin_complementaire:" "Assertions effectuées avec succès."