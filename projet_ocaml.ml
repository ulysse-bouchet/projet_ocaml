(*
@file projet_ocaml.ml
@author Ulysse Bouchet
@version 0.5.1
@data 13/05/2021
*)

(* Imports *)
open Format

(*************************)
(*** SÉQUENCES ADN/ARN ***)
(*************************)

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
let contenu_gc (brin: brin): float =
    if List.length brin = 0 then 
        failwith "[contenu_gc] erreur : brin vide"
    else
        List.fold_left (fun acc nucl -> acc +. if nucl = G || nucl = C then 1. else 0.) 0. brin
        /. 
        (float) (List.length brin)
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
    let rec brin_vers_chaine_rec (b: brin) (chaine: chaine) (start: bool): chaine =
        match b with
            | []            -> chaine
            | _::[]         -> failwith "[brin_vers_chaine] erreur : brin invalide"
            | _::_::[]      -> failwith "[brin_vers_chaine 188] erreur : brin invalide"
            | n1::n2::n3::q -> 
                let acide = 
                    codon_vers_acide n1 n2 n3
                in 
                    match (q, start, acide) with
                    | (_, true, STOP)       -> chaine (* On arrive au bout du brin *)
                    | ([], _, _)            -> failwith "[brin_vers_chaine] erreur : brin invalide" (* On arrive au bout du brin sans STOP => brin invalide *)
                    | (_, false, START)     -> brin_vers_chaine_rec q [] true (* On commence la chaîne *)
                    | (_, true, START)      -> failwith "[brin_vers_chaine] erreur : brin invalide" (* On a déjà commencé la chaîne => brin invalide *)
                    | (_, false, _)         -> failwith "[brin_vers_chaine] erreur : brin invalide" (* On a un acide sans avoir commencé de chaîne => brin invalide *)
                    | (_, true, _)          -> brin_vers_chaine_rec q (chaine@[acide]) true (* On passe au codon suivant *)
    in
        brin_vers_chaine_rec b [] false
    ;;

(* Assertions *)
let () = assert (brin_vers_chaine [T;A;C;G;G;C;T;A;G;A;T;T;T;A;C;G;C;T;A;A;T;A;T;C] = [Pro;Ile]);;
let () = assert (brin_vers_chaine [] = []);;
(* brin_vers_chaine [T;A;C;T;A;C] provoquerait une erreur car on y trouve deux acides START d'affilée *)
(* brin_vers_chaine [T;A;C;G;G;A;T;C] provoquerait une erreur car la taille du brin n'est pas un multiple de 3 (donc indivisible en codons) *)
(* brin_vers_chaine [G;G;C;A;T;T] provoquerait une erreur car le brin ne commence pas par START*)
(* brin_vers_chaine [T;A;C;G;G;C] provoquerait une erreur car le brin ne commence pas par START*)
let () = printf "%-30s %s\n" "brin_vers_chaine:" "Assertions effectuées avec succès.";;

(*** Question 6 ***)

(*
@function brin_vers_chaine
@desc Décode toutes les chaînes d'un brin d'ADN donné.
@param b Un brin d'ADN
@returns Toutes les chaînes d'acides aminés du brin d'ADN 
@throws Une exception lorsque le brin est invalide
*)
let brin_vers_chaines (b: brin): chaine list = 
    let rec brin_vers_chaines_rec (b: brin) (chaine: chaine) (chaines: chaine list) (start: bool): chaine list =
        match b with
            | []            -> chaines
            | _::[]         -> failwith "[brin_vers_chaines] erreur : brin invalide"
            | _::_::[]      -> failwith "[brin_vers_chaines] erreur : brin invalide"
            | n1::n2::n3::q -> 
                let acide = 
                    codon_vers_acide n1 n2 n3
                in 
                    match (q, start, acide) with
                    | ([], true, STOP)      -> chaines@[chaine] (* On arrive au bout du brin *)
                    | ([], _, _)            -> failwith "[brin_vers_chaines] erreur : brin invalide" (* On arrive au bout du brin sans STOP => brin invalide *)
                    | (_, false, START)     -> brin_vers_chaines_rec q [] chaines true (* On commence une nouvelle chaîne *)
                    | (_, true, START)      -> failwith "[brin_vers_chaines] erreur : brin invalide" (* On a déjà commencé une chaîne => brin invalide *)
                    | (_, true, STOP)       -> brin_vers_chaines_rec q [] (chaines@[chaine]) false (* On passe à la chaîne suivante *)
                    | (_, false, _)         -> failwith "[brin_vers_chaines] erreur : brin invalide" (* On a un acide sans avoir commencé de chaîne => brin invalide *)
                    | (_, true, _)          -> brin_vers_chaines_rec q (chaine@[acide]) chaines true (* On passe au codon suivant *)
    in
        brin_vers_chaines_rec b [] [] false
    ;;

let () = assert (brin_vers_chaines [T;A;C;G;G;C;T;A;G;A;T;T;T;A;C;G;C;T;A;A;T;A;T;C] = [[Pro;Ile]; [Arg;Leu]]);;
let () = assert (brin_vers_chaines [] = []);;
(* brin_vers_chaines [T;A;C;G;G;C;T;A;G;A;T;T;T;A;C;G;C;T;A;A;T];; provoquerait une erreur car ne termine pas par STOP*)
(* brin_vers_chaines [T;A;C;T;A;C;T;A;G;A;T;T;T;A;C;G;C;T;A;A;T;A;T;C];; provoquerait une erreur car possède deux START sans STOP au milieu*)
(* brin_vers_chaines [G;G;C;T;A;G;A;T;T;T;A;C;G;C;T;A;A;T;A;T;C];; provoquerait une erreur car ne commence pas par START*)
(* brin_vers_chaines [T;A;C;G;G;C;T;A;G;A;T;T;G;C;T;A;A;T;A;T;C];; provoquerait une erreur car une chaine commence sans START*)
let () = printf "%-30s %s\n" "brin_vers_chaines:" "Assertions effectuées avec succès.";;

(******************************)
(*** ARBRES PHYLOGÉNÉTIQUES ***)
(******************************)

(* Arbres phylogénétiques - Type de base *)
(* @author Stefania Dumbrava *)
type arbre_phylo = 
    | Lf of brin 
    | Br of arbre_phylo * brin * int * arbre_phylo;;

(*** Question 1 ***)

(*
@function brin_vers_string
@desc Donne la représentation en chaîne de caractères d'un brin donné (ex: [A;T;C;G] --> ATCG)
@param Un brin d'ADN
@returns Sa représentation en chaîne de caractères
*)
let brin_vers_string (brin: brin) =
    List.fold_left (fun str nucl -> str ^ match nucl with |A->"A"|T->"T"|C->"C"|G->"G") "" brin
    ;;

(*
@function arbre_phylo_vers_string
@desc Donne la représentation en chaîne de caractères d'un arbre phylogénétique donné
@param Un arbre phylogénétique
@returns Sa représentation en chaîne de caractères
*)
let rec arbre_phylo_vers_string (arbre: arbre_phylo): string =
    match arbre with 
    | Lf (brin) ->
        "(" ^ brin_vers_string brin ^ ")"
    | Br (arbre_gauche, brin, malus, arbre_droit) ->
        "{" ^ arbre_phylo_vers_string arbre_gauche ^ "}" ^
        " <-- (" ^ brin_vers_string brin ^ " " ^ string_of_int malus ^ ") --> " ^
        "{" ^ arbre_phylo_vers_string arbre_droit ^ "}"
    ;;

(* Définition de quelques arbres *)

let arbre_1 = Br (Br (Lf ([G;C;A;T]), [A;C;A;T], 3, Lf ([T;C;G;T])), [A;A;A;A], 8, Br (Lf ([T;A;G;A]), [A;A;G;A], 2, Lf ([G;A;G;A])));;
let arbre_2 = Br (Br (Lf ([G;A;A;T]), [G;C;T;T], 5, Lf ([C;A;G;T])), [A;A;T;A], 12, Br (Lf ([T;A;G;A]), [A;A;G;A], 3, Lf ([G;T;G;A])));;
let arbre_3 = Lf ([]);;
let arbre_4 = Br (Br (Lf ([G;A;G;T]), [G;C;T;T], 5, Lf ([C;A;G;T])), [A;A;T;A], 12, Br (Lf ([T;A;G;A]), [A;A;G;A], 3, Lf ([G;T;G;A])));;

(* Assertions *)

let () = assert (arbre_phylo_vers_string arbre_1 = 
                    "{{(GCAT)} <-- (ACAT 3) --> {(TCGT)}} <-- (AAAA 8) --> {{(TAGA)} <-- (AAGA 2) --> {(GAGA)}}");;
let () = assert (arbre_phylo_vers_string arbre_2 = 
                    "{{(GAAT)} <-- (GCTT 5) --> {(CAGT)}} <-- (AATA 12) --> {{(TAGA)} <-- (AAGA 3) --> {(GTGA)}}");;
let () = assert (arbre_phylo_vers_string arbre_3 = 
                    "()");;
let () = printf "%-30s %s\n" "arbre_phylo_vers_string:" "Assertions effectuées avec succès.";;

(*** Question 2 ***)

(*
@function similarite_arbre
@desc Calcule la similarité de deux arbres phylogénétiques en comparant leurs éléments deux à deux.
@param arbre_1 Un arbre phylogénétique 
@param arbre_2 Un autre arbre phylogénétique
@returns La somme de leur similarités procentuelles
*)
let rec similarite_arbre (arbre_1: arbre_phylo) (arbre_2: arbre_phylo): float =
    match (arbre_1, arbre_2) with
    | Lf (_), Br (_, _, _, _) ->
        failwith "[similarite_arbre] erreur : tailles différentes"
    | Br (_, _, _, _), Lf (_) ->
        failwith "[similarite_arbre] erreur : tailles différentes"
    | Lf (brin_1), Lf (brin_2) ->
        similarite brin_1 brin_2
    | Br (arbre_gauche_1, brin_1, malus_1, arbre_droit_1), Br (arbre_gauche_2, brin_2, malus_2, arbre_droit_2) ->
        (similarite brin_1 brin_2) +. (similarite_arbre arbre_gauche_1 arbre_gauche_2) +. (similarite_arbre arbre_droit_1 arbre_droit_2)
    ;;

(*
@function similaire
@desc Calcule, pour une liste d'arbres et un arbre donné, l'arbre le plus similaire avec ce dernier.
@param arbre Un arbre phylogénétique
@param arbres Une liste d'arbres phylogénétiques
@returns L'arbre le plus similaire
*)
let similaire (arbre: arbre_phylo) (arbres: arbre_phylo list): arbre_phylo =
    let rec aux restants plus_grand plus_grand_simi =
        match restants with
        | []    -> plus_grand
        | t::q  -> let simi = similarite_arbre arbre t in if simi < plus_grand_simi then aux q plus_grand plus_grand_simi else aux q t simi
    in
        match arbres with
        | []    -> failwith "[similaire] erreur : liste vite"
        | t::q  -> aux q t (similarite_arbre arbre t)
    ;;

(* Assertions *)

let () = assert (similaire arbre_1 [arbre_4; arbre_2] = arbre_2);;
let () = assert (similaire arbre_1 [arbre_1; arbre_2; arbre_4] = arbre_1);;
(* similaire arbre_1 [] provoquerait une exception puisque la liste est vide *)
let () = printf "%-30s %s\n" "similaire:" "Assertions effectuées avec succès.";;

(*** Question 3 ***)

(*
@function get_root
@desc Extrait le brin de la racine d'un arbre phylogénétique.
@param arbre Un arbre phylogénétique
@returns Le brin de sa racine
*)
let get_root (arbre: arbre_phylo): brin =
    match arbre with
    | Lf (brin)             -> brin
    | Br (_, brin, _, _)    -> brin
    ;;

(* Assertions *)

let () = assert (get_root arbre_1 = [A;A;A;A]);;
let () = assert (get_root arbre_2 = [A;A;T;A]);;
let () = assert (get_root arbre_3 = []);;
let () = printf "%-30s %s\n" "get_root:" "Assertions effectuées avec succès.";;

(*
@function get_malus
@desc Extrait le malus de la racine d'un arbre phylogénétique.
@param arbre Un arbre phylogénétique
@returns Le malus de sa racine
*)
let get_malus (arbre: arbre_phylo): int =
    match arbre with
    | Lf (_)                    -> 0
    | Br (_, _, malus, _)       -> malus
    ;;

(* Assertions *)

let () = assert (get_malus arbre_1 = 8);;
let () = assert (get_malus arbre_2 = 12);;
let () = assert (get_malus arbre_3 = 0);;
let () = printf "%-30s %s\n" "get_malus:" "Assertions effectuées avec succès.";;

(*
@function br
@desc Construit un arbre phylogénétique à partir d'un sous-arbre gauche, d'un brin et d'un sous-arbre droit.
@param arbre_gauche Un arbre phylogénétique
@param brin Un brin d'ADN
@param arbre_droit Un autre arbre phylogénétique
@returns Un arbre phylogénétique basé sur les paramètres
*)
let br (arbre_gauche: arbre_phylo) (brin: brin) (arbre_droit: arbre_phylo): arbre_phylo =
    let malus =
        (distance brin (get_root arbre_gauche) + get_malus arbre_gauche) + (distance brin (get_root arbre_droit) + get_malus arbre_droit)
    in
        Br (arbre_gauche, brin, malus, arbre_droit)
    ;;

(* Assertions *)

let br_res = Br (Br (Lf ([A;A;T;T]), [T;T;T;T], 4, Lf ([T;T;A;A])), [A;A;A;A], 14, Br (Lf ([A;A;T;T]), [T;A;A;T], 4, Lf ([T;T;A;A])));;
let br_brin = [A;A;A;A];;
let br_gauche = Br (Lf ([A;A;T;T]), [T;T;T;T], 4, Lf ([T;T;A;A]));;
let br_droit = Br (Lf ([A;A;T;T]), [T;A;A;T], 4, Lf ([T;T;A;A]));;

let () = assert (br br_gauche br_brin br_droit = br_res);;
let () = printf "%-30s %s\n" "br:" "Assertions effectuées avec succès.";;

(*** Question 4 ***)

(*
@function gen_phylo
@desc Génère tous les arbres phylogénétiques possibles à partir des trois brins donnés.
@param b1 Un brin d'ADN
@param b2 Un autre brin d'ADN
@param b3 Un autre brin d'ADN
@returns La liste des arbres phylogénétiques générés
*)
let gen_phylo (b1: brin) (b2: brin) (b3: brin): arbre_phylo list =
    [
        br (Lf (b1)) b2 (Lf (b3));
        br (Lf (b1)) b3 (Lf (b2));
        br (Lf (b2)) b1 (Lf (b3));
        br (Lf (b2)) b3 (Lf (b1));
        br (Lf (b3)) b1 (Lf (b2));
        br (Lf (b3)) b2 (Lf (b1))
    ]
    ;;

(* Assertions *)

let b1 = [A;A;A;A];;
let b2 = [T;T;T;T];;
let b3 = [C;C;C;C];;
let a1 = Br (Lf (b1), b2, distance b2 b1 + distance b2 b3, Lf (b3));;
let a2 = Br (Lf (b1), b3, distance b3 b1 + distance b3 b2, Lf (b2));;
let a3 = Br (Lf (b2), b1, distance b1 b2 + distance b1 b3, Lf (b3));;
let a4 = Br (Lf (b2), b3, distance b3 b2 + distance b3 b1, Lf (b1));;
let a5 = Br (Lf (b3), b1, distance b1 b3 + distance b1 b2, Lf (b2));;
let a6 = Br (Lf (b3), b2, distance b2 b3 + distance b2 b1, Lf (b1));;
let a_res = [a1; a2; a3; a4; a5; a6];;

let () = assert (gen_phylo b1 b2 b3 = a_res);;
let () = printf "%-30s %s\n" "gen_phylo:" "Assertions effectuées avec succès.";;

(*** Question 5 ***)

(*
@function min_malus
@desc Pour une liste donnée, renvoie l'arbre phylogénétique ayant le malus minimal.
@param arbres Une liste d'arbres phylogénétiques
@returns L'arbre ayant le malus minimal
*)
let min_malus (arbres: arbre_phylo list): arbre_phylo =
    let rec min_malus_rec (arbres: arbre_phylo list) (min_malus: int) (min_arbre: arbre_phylo): arbre_phylo =
        match arbres with
        | []    -> min_arbre
        | t::q  -> 
            let cur_malus = 
                get_malus t 
            in 
                min_malus_rec q (if cur_malus < min_malus then cur_malus else min_malus) (if cur_malus < min_malus then t else min_arbre)
    in
        match arbres with
        | []    -> failwith "[min_malus] erreur : liste vide" 
        | t::q  -> min_malus_rec q (get_malus t) t
    ;;

(* Assertions *)

let () = assert (min_malus [arbre_1; arbre_2] = arbre_1);;
let () = printf "%-30s %s\n" "min_malus:" "Assertions effectuées avec succès.";;
