(*
@file projet_ocaml.ml
@author Ulysse Bouchet
@version 0.5.1
@data 13/05/2021
*)

(* Imports *)
open Format

(* Séquences génétiques - Types de bases *)
(* @author : Stefania Dumbrava *)
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
        (float) (List.length b)
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

(*** Question 4 ***)

(* 
@function similarite
@desc Calcule la similarité procentuelle entre deux brins de même longueur.
@param b1 Un brin d'ADN
@param b2 Un autre brin d'ADN
@returns La similarité entre les brins
@see distance
*)
let similarite (b1: brin) (b2: brin): float =
    if b1 = [] && b2 = [] then
        1.
    else 
        let length = 
            (float) (List.length b1)
        in 
            1. -. ((float) (distance b1 b2) /. length)
    ;;

(* Assertions *)
let () = assert (similarite [C;G;A;T] [T;A;G;T] = 0.25);;
let () = assert (similarite [A;G;C;T] [T;A;A;G] = 0.);;
let () = assert (similarite [A;G;C;T] [A;G;C;T] = 1.);;
let () = assert (similarite [] [] = 1.);;
(* similarite [A] [A;A] provoquerait une erreur car les listes sont de longueur différentes*)
let () = printf "%-30s %s\n" "similarite:" "Assertions effectuées avec succès.";;

(*** Question 5 ***)

(* Séquences génétiques - Type avancé : acide aminé *)
(* @author : Stefania Dumbrava *)
type acide =    Ala | Arg | Asn | Asp | Cys | Glu | Gln | Gly | His | Ile | 
                Leu | Lys | Phe | Pro | Ser | Thr | Trp | Tyr | Val | 
                START | STOP;;

(* Séquences génétiques - Type avancé : chaine *)
type chaine = acide list;;

(*
@function codon_vers_acide
@desc Pour un codon donné, renvoie l'acide aminé qui lui correspond (ou START/STOP).
@param n1 Le premier nucléotide du codon
@param n2 Le deuxième nucléotide du codon
@param n3 Le troisième nucléotide du codon
@returns L'acide aminé correspondant au codon donné.
@author : Stefania Dumbrava
*)
let codon_vers_acide n1 n2 n3 = 
    match n1, n2, n3 with
    | (A,A,A) -> Phe  | (A,A,G) -> Phe  | (A,A,T) -> Leu  | (A,A,C) -> Leu 
    | (G,A,A) -> Leu  | (G,A,G) -> Leu  | (G,A,T) -> Leu  | (G,A,C) -> Leu
    | (T,A,A) -> Ile  | (T,A,G) -> Ile  | (T,A,T) -> Ile  | (T,A,C) -> START
    | (C,A,A) -> Val  | (C,A,G) -> Val  | (C,A,T) -> Val  | (C,A,C) -> Val 
    | (A,G,A) -> Ser  | (A,G,G) -> Ser  | (A,G,T) -> Ser  | (A,G,C) -> Ser
    | (G,G,A) -> Pro  | (G,G,G) -> Pro  | (G,G,T) -> Pro  | (G,G,C) -> Pro
    | (T,G,A) -> Thr  | (T,G,G) -> Thr  | (T,G,T) -> Thr  | (T,G,C) -> Thr
    | (C,G,A) -> Ala  | (C,G,G) -> Ala  | (C,G,T) -> Ala  | (C,G,C) -> Ala
    | (A,T,A) -> Tyr  | (A,T,G) -> Tyr  | (A,T,T) -> STOP | (A,T,C) -> STOP
    | (G,T,A) -> His  | (G,T,G) -> His  | (G,T,T) -> Gln  | (G,T,C) -> Gln
    | (T,T,A) -> Asn  | (T,T,G) -> Asn  | (T,T,T) -> Lys  | (T,T,C) -> Lys
    | (C,T,A) -> Asp  | (C,T,G) -> Asp  | (C,T,T) -> Glu  | (C,T,C) -> Glu
    | (A,C,A) -> Cys  | (A,C,G) -> Cys  | (A,C,T) -> STOP | (A,C,C) -> Trp
    | (G,C,A) -> Arg  | (G,C,G) -> Arg  | (G,C,T) -> Arg  | (G,C,C) -> Arg
    | (T,C,A) -> Ser  | (T,C,G) -> Ser  | (T,C,T) -> Arg  | (T,C,C) -> Arg
    | (C,C,A) -> Gly  | (C,C,G) -> Gly  | (C,C,T) -> Gly  | (C,C,C) -> Gly
    ;;

(*
@function brin_vers_chaine
@desc Décode les codons de la première chaîne d'acide aminés d'un brin d'ADN donné.
@param b Un brin d'ADN
@returns La première chaîne d'acide aminés du brin
@throws Une exception lorsque le brin est invalide
*)
let brin_vers_chaine (b: brin): chaine =
    let rec brin_vers_chaine_rec (b: brin) (chaine: chaine): chaine =
        match b with
            | []            -> failwith "[brin_vers_chaine] erreur : brin invalide" (* Si ce cas est atteint, alors il n'y a pas eu de STOP, le brin est donc invalide *)
            | _::[]         -> failwith "[brin_vers_chaine] erreur : brin invalide"
            | _::_::[]      -> failwith "[brin_vers_chaine] erreur : brin invalide"
            | n1::n2::n3::q -> 
                let acide = 
                    codon_vers_acide n1 n2 n3
                in 
                    if acide = STOP then
                        chaine
                    else if acide = START then
                        failwith "[brin_vers_chaine] erreur : brin invalide"
                    else
                        brin_vers_chaine_rec q (chaine@[acide])
    in
        match b with
        | []            -> []
        | _::[]         -> failwith "[brin_vers_chaine] erreur : brin invalide"
        | _::_::[]      -> failwith "[brin_vers_chaine] erreur : brin invalide"
        | n1::n2::n3::q -> 
            let acide = 
                codon_vers_acide n1 n2 n3
            in 
                if acide = START then
                    brin_vers_chaine_rec q []
                else 
                    failwith "[brin_vers_chaine] erreur : brin invalide"
    ;;

(* Assertions *)
let () = assert (brin_vers_chaine [T;A;C;G;G;C;T;A;G;A;T;T;T;A;C;G;C;T;A;A;T;A;T;C] = [Pro;Ile]);;
let () = assert (brin_vers_chaine [] = []);;
(* brin_vers_chaine [T;A;C;T;A;C] provoquerait une erreur car on y trouve deux acides START d'affilée *)
(* brin_vers_chaine [T;A;C;G;G;A;T;C] provoquerait une erreur car la taille du brin n'est pas un multiple de 3 (donc indivisible en codons) *)
(* brin_vers_chaine [G;G;C;A;T;T] provoquerait une erreur car le brin ne commence pas par START*)
(* brin_vers_chaine [T;A;C;G;G;C] provoquerait une erreur car le brin ne commence pas par START*)
let () = printf "%-30s %s\n" "brin_vers_chaine:" "Assertions effectuées avec succès.";;