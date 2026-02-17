;; Initial ns statement copied from
;; https://generateme.github.io/fastmath/clay/vector_matrix.html#matrices
(ns scratch
  (:require ;[clojure.math :as m]
            ;[criterium.core :as crit]
            [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]
            [fastmath.random :as frand]
            [scicloj.kindly.v4.kind :as kind]
            [scicloj.tableplot.v1.plotly :as plotly]
            [tablecloth.api :as tc]
            [drifthappens.wrightfisher :as wf]
            [utils.plotly :as uplot]
            [utils.fastmats :as mats]
            [utils.misc :as umisc])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))

;; Example matrix from https://generateme.github.io/fastmath/clay/vector_matrix.html
(def M3x3 (Mat3x3. 1 2 3 -4 5 6 9 -8 7))

(def m3m (fmat/mulm M3x3 M3x3))

;; Generates a structure containing a map in which many components
;; are not delayed and not yet forced.  fmat/decomposition-component 
;; forces and replaces the :not-delivered value with the actual value.
(def m3m-edecomp (fmat/eigen-decomposition m3m))

(def v (fmat/decomposition-component m3m-edecomp :V))
(def d (fmat/decomposition-component m3m-edecomp :D))
(def vi (fmat/inverse v))

;; mulm takes only two matrix args
(reduce fmat/mulm [v d vi]) ; this product should be equal to m3m (up to float-slop)

v
;; A sequence of the columns of v:
(fmat/eigenvectors m3m)
(fmat/eigenvectors m3m true)

d
(fmat/eigenvalues m3m) ; seq of pairs of [real, imag] components of eigenvals

;; Apparently these are not real and imaginary evals, but rather
;; real and imaginary components of eigenvalues:
(fmat/decomposition-component m3m-edecomp :real-eigenvalues) ; ? diag of d
(fmat/decomposition-component m3m-edecomp :imag-eigenvalues) ; ? sub/super pairs of d
(fmat/decomposition-component m3m-edecomp :eigenvectors) ; ? sub/super pairs of d

;; Seems to be same as fmat/eigenvalues but as a matrix rather than
;; a clojure sequence.  Also seems to be same as :D from the
;; eigen-decomposition.
(fmat/eigenvalues-matrix m3m)

;; Seems to be same as :V in the eigen-decomposition:
(fmat/eigenvectors m3m)

;; --------------------------

;; From predreprod1.clj:

;; ---
;; ### Configuration and setup

(def fit-B 1.0)
(def big-fit-A 1.01)  ; large size has sel
(def small-fit-A 1.0) ; small size is pure drift
(def interval 16) ; interval between generations to display
(def half-interval (/ interval 2))

(def increments (iterate (partial + 16) 0))
(def generations increments)
(def half-generations (map (fn [n] (/ n 2)) increments))

(def num-gens 6)

;; These next two values should even numbers so that when divided,
;; we'll get integers:
(def small-N 6)
(def big-N 14)

(def half-small-N (/ small-N 2))
(def half-big-N (/ big-N 2))

;; ---
;; ### Small populations

(def small-pop-init 
  "A vector of initial probabilities of frequencies for a small population,
  typically with 1 as one element and zeros elswhere."
  (mats/mkvec (concat (repeat half-small-N 0.0) [1.0] (repeat half-small-N 0.0))))

(def small-tran-mat 
  "A single transition matrix for a small population."
  (wf/right-mult-tran-mat small-fit-A fit-B (dec (count small-pop-init))))

;; NOTE: Consider replacing choose-mat-powers-separately with choose-mat-powers-sequentially if the exponents are closely spaced; it might be more efficient:
(def small-tran-mats 
  "A sequence of M-to-M transition matrices, each of which is small-tran-mat raised to a power."
  (doall (mats/choose-mat-powers-separately small-tran-mat (take num-gens generations))))

(def small-prob-states 
  "States resulting from applying a product transition matrix to an initial state."
  (mats/make-prob-states small-tran-mats small-pop-init))

;; Plots made from the preceding sequence of states.
(def small-plots (mapv uplot/plot-both small-prob-states))

;small-plots  ; display the plots

;; ---
;; ### Large populations

(def big-pop-init 
  "A vector of initial probabilities of frequencies for a big population,
  typically with 1 as one element and zeros elswhere."
  (mats/mkvec (concat (repeat half-big-N 0.0) [1.0] (repeat half-big-N 0.0))))

(def big-tran-mat 
  "A single transition matrix for a big population."
  (wf/right-mult-tran-mat big-fit-A fit-B (dec (count big-pop-init)))) ; use fit-B for fit-A to make them equal

;; NOTE: Consider replacing choose-mat-powers-separately with choose-mat-powers-sequentially if the exponents are closely spaced; it might be more efficient:
(def big-tran-mats
  "A sequence of N-to-N transition matrices, each of which is big-tran-mat raised to a power."
  (doall (mats/choose-mat-powers-separately big-tran-mat (take num-gens generations))))

(def big-prob-states 
  "States resulting from applying a product transition matrix to an initial state."
  (mats/make-prob-states big-tran-mats big-pop-init))

;; Plots made from the preceding sequence of states.
(def big-plots (mapv uplot/plot-lines big-prob-states))

;big-plots  ; display the plots

;; ---
;; ### Alternating small/large populations

(def N (dec (count big-pop-init)))
(def M (dec (count small-pop-init)))

;; Shrinking and expanding matrices.
;;
;; Note M and N are swapped in the next two
(def predat-tran-mat
  "A transition matrix from a population of size N to a smaller population
  of size M resulting from predation."
  (wf/right-mult-tran-mat small-fit-A fit-B N M)) ; 

(def reprod-tran-mat
  "A transition matrix from a population of size M to a larger population
  of size N due to reproduction."
  (wf/right-mult-tran-mat big-fit-A fit-B M N)) ; use fit-B for fit-A to make them equal

;; Combine into one square matrix:
(def pred-reprod-mat 
  "A transition matrix from a populatiion of size N to a population of the
  same size. The matrix is the product of predat-tran-mat and reprod-drift
  mat.  It represents predation followed by reproductive growth."
  (fmat/mulm reprod-tran-mat predat-tran-mat))

;; The problem below with eigenvectors is this one:
;; https://github.com/generateme/fastmath/issues/42

(comment
  ;; experiments with eigenvecs/vals
  pred-reprod-mat
  (type pred-reprod-mat)
  (fmat/shape pred-reprod-mat)

  ;; I get "illegal state: convergence failed" errors, even for 100x100 matrix
  (fmat/eigenvalues pred-reprod-mat) ;; works with 5x5
  (fmat/eigenvectors pred-reprod-mat) ; fails.  Tried various conversions
  (def predreproddecomp (fmat/eigen-decomposition pred-reprod-mat))
  (fmat/mat->array 
    (fmat/decomposition-component predreproddecomp :D)
  )
  (fmat/mat->array 
    (fmat/decomposition-component predreproddecomp :V)
  )
  (type (fmat/decomposition-component predreproddecomp :V))
  (fmat/decomposition-component predreproddecomp :real-eigenvalues)
  (fmat/decomposition-component predreproddecomp :imag-eigenvalues)
  (fmat/decomposition-component predreproddecomp :eigenvectors)
  (map fvec/vec->seq (fmat/decomposition-component predreproddecomp :eigenvectors))

  ;; ---

(def m (fmat/rows->RealMatrix
           [[1.0, 0.6637495573585641, 0.4272874559265132, 0.2654410598289928,
             0.15811598711862837, 0.08956507176600371, 0.0477138732755275,
             0.023542733926677765, 0.010525383640630919, 0.004124092472272302,
             0.001341370727588796, 3.2821670644383045E-4, 4.8912070734631946E-5,
             2.3648379316104523E-6, 0.0],
            [0.0, 0.06498727081139731, 0.08918591901162025, 0.08899340480302088,
             0.076094979818456, 0.05837576249647409, 0.0407409763701028,
             0.02584435126923635, 0.014724687436623222, 0.00735058255745415,
             0.0030733217025500346, 9.87930185150087E-4, 2.023893313002936E-4,
             1.5015163841875311E-5, 0.0],
            [0.0, 0.08794333963174404, 0.12477929280788937, 0.1292039760036972,
             0.11514586577659561, 0.0925533100420665, 0.06811876562294183,
             0.04593976207210663, 0.028116591768610897, 0.01528672628718937,
             0.007095959041189156, 0.002606274198905312, 6.404418733248367E-4,
             6.333958527836005E-5, 0.0],
            [0.0, 0.07639058338382698, 0.1166468933399246, 0.13017778015811657,
             0.1252800345401505, 0.1090186925379813, 0.08715620175421784,
             0.0641304488047315, 0.04308164804342768, 0.025928091549180268,
             0.01349076037492816, 0.00566679705893444, 0.0016518393982085998,
             2.1121548409098616E-4, 0.0],
            [0.0, 0.04997940948378009, 0.08765529258270305, 0.11042273537781772,
             0.11851032053067201, 0.11406117428824489, 0.10033772699230568,
             0.08105018709319961, 0.059808728668060064, 0.039699392443448216,
             0.022983700322418473, 0.010921983416010885, 0.0037204235791706373,
             6.018084510938057E-4, 0.0],
            [0.0, 0.028094043799514377, 0.06029808517253643, 0.0878429223155858,
             0.10587050564799562, 0.11258706062787878, 0.1084844915575997,
             0.09563006637875067, 0.07702438245663279, 0.05602761335424251,
             0.03585403659576232, 0.019134842419556316, 0.00754922352067117,
             0.0015237457828414425, 0.0],
            [0.0, 0.014930908516603777, 0.039906515305181745, 0.06680593112294475,
             0.08968775592165211, 0.1047162211065144, 0.1100204313409928,
             0.10548670589749444, 0.09248401955396415, 0.07352254303637273,
             0.05184528300710122, 0.030952821599221296, 0.014061155496672447,
             0.0034926345603350627, 0.0],
            [0.0, 0.0076605496089470655, 0.025054877280984794,
             0.04780303995739107, 0.07130684615121308, 0.09137076092832726,
             0.10468184940987312, 0.1091487587960605, 0.10409973891135059,
             0.09033970127101076, 0.07006631666904285, 0.046645151287484844,
             0.02424384132708681, 0.007325306159360033, 0.0],
            [0.0, 0.003684123852951813, 0.014642409058654993, 0.03194110059528582,
             0.05311005047358716, 0.07483765144474376, 0.0935901904512271,
             0.10615562182060892, 0.11011176020234424, 0.10421889324752348,
             0.08873681403159372, 0.06566627322004881, 0.03891485097708908,
             0.014387248617250003, 0.0],
            [0.0, 0.0016204846972839295, 0.007915543801113114,
             0.019872979432674524, 0.036958226329897544, 0.057381953876503265,
             0.07842778050521788, 0.09684010047514993, 0.1092720230233315,
             0.1127934238904237, 0.10545910922058574, 0.08693709183550805,
             0.05919697988260963, 0.027258477857399113, 0.0],
            [0.0, 6.445970134384694E-4, 0.003925367140393781,
             0.011410646075883891, 0.023828693952606724, 0.04089269979843884,
             0.061249484431984524, 0.08255190179617528, 0.10165493872991943,
             0.11493551317780272, 0.11873597083783832, 0.10993128024726814,
             0.08661992630641394, 0.04893850224057894, 0.0],
            [0.0, 2.2764395941963882E-4, 0.0017523145006021325,
             0.0059503369768318356, 0.014054998321100099, 0.02683518385157636,
             0.04432874325163541, 0.06563058572444108, 0.08873550432208584,
             0.11043572944928717, 0.12627321154163942, 0.13054663291842283,
             0.11637314880996794, 0.07580485755957672, 0.0],
            [0.0, 6.861135452365823E-5, 6.82135898864356E-4,
             0.0027467130278743417, 0.007419728778561235, 0.015883384642737815,
             0.029058804158331366, 0.04726067727935268, 0.0697924425245244,
             0.09448200690456945, 0.11715800362815809, 0.1310665875865158,
             0.12622876861668958, 0.08873828254347459, 0.0],
            [0.0, 1.631523517687666E-5, 2.1610012014743415E-4,
             0.0010444190384004516, 0.0032274382603528435, 0.0076817619686832235,
             0.015330927358408552, 0.026829058811703943, 0.04219168114746543,
             0.06033369194561612, 0.07851449294615492, 0.09169028052294907,
             0.09177349523226917, 0.06679943043606759, 0.0],
            [0.0, 2.5612928278137376E-6, 5.1798052871113494E-5,
             3.4295528548221993E-4, 0.001388568378531129, 0.004239310623825501,
             0.010759753519633706, 0.023959039854311047, 0.048376469571028806,
             0.09052199841360695, 0.15937164935344883, 0.26691783679758035,
             0.42877460357779124, 0.6648377707208799, 1.0]]))

(= pred-reprod-2d m)

  ;; What happens if I convert to RealMat first?
  (def pred-reprod-2d (fmat/array2d->RealMatrix pred-reprod-mat))
  ;(def pred-reprod-2d (fmat/array2d->RealMatrix big-tran-mat)) ;; FIXME TEST
  (type pred-reprod-2d)
  (fmat/shape pred-reprod-2d)
  (fmat/mat->array2d pred-reprod-2d)

  (def ra2dcomp (fmat/eigen-decomposition pred-reprod-2d))
  (def evals (fmat/decomposition-component ra2dcomp :real-eigenvalues)) ; no imaginary parts in this case 
  (def evecs (fmat/decomposition-component ra2dcomp :eigenvectors))
  (fmat/eigenvectors pred-reprod-2d) ; fails.  Tried various conversions
  (fmat/decomposition-component ra2dcomp :D)
  (fmat/decomposition-component ra2dcomp :V)
  (type (fmat/decomposition-component ra2dcomp :V))
  (type evals) ; evals is Java vector not a Clojure vector
  (fmat/decomposition-component ra2dcomp :imag-eigenvalues)
  (map fvec/vec->seq evecs)

  ;; test the eigenvector/value relationship for each pair:
  (map #(fvec/delta-eq (fmat/mulv pred-reprod-2d %1) ; matrix-vector mult
                       (fvec/mult %1 %2)) ; vector-scalar mult
       evecs
       (vec evals))


  (map (fn [evec] fmat/mulv pred-reprod-2d evec) evecs)
  (map (fn [evec] fmat/mulv pred-reprod-2d evec) evecs)



)

(comment
;; We use half-generations here because each step involves two sampling processes.  
;; So each generation is analogous to two generations in the small and big models.
;; NOTE: Consider replacing choose-mat-powers-separately with choose-mat-powers-sequentially if the exponents are closely spaced; it might be more efficient:
(def pred-reprod-tran-mats
  "A sequence of N-to-N transition matrices, each of which is pred-reprod-mat raised to a power."
  (doall (mats/choose-mat-powers-separately pred-reprod-mat (take num-gens half-generations))))

(def pred-reprod-prob-states
  "States resulting from applying a product transition matrix to an initial state."
  (mats/make-prob-states pred-reprod-tran-mats big-pop-init))

;; Plots made from the preceding sequence of states.
(def pred-reprod-plots (mapv uplot/plot-both pred-reprod-prob-states))
)

(clojure.repl/dir fastmath.matrix)

;pred-reprod-plots  ; display the plots


;; ---
; ### Plots from all three populations:

(comment
(def small-big-combo-plots 
  "States resulting from alternating plots from other sequences of plots."
  (interleave small-plots big-plots pred-reprod-plots))
)

;small-big-combo-plots   ; display the plots


;; ----------------------------------------

;; double-double, toil and trouble
(def witch (mats/mkmat [[1 2 3 4]
                        [5 6 7 8]
                        [9 -1 -2 -3]
                        [-4 -5 -6 -7]]))

;; Note this takes four seqs, not a seq of seqs:
(def truck (fmat/rows->mat4x4 [1 2 3 4]
                              [5 6 7 8]
                              [9 -1 -2 -3]
                              [-4 -5 -6 -7]))

(def witchdecomp (fmat/eigen-decomposition witch))
(def truckdecomp (fmat/eigen-decomposition truck))

(def tv (fmat/decomposition-component witchdecomp :V))
(def wv (fmat/decomposition-component truckdecomp :V))
(def td (fmat/decomposition-component witchdecomp :D))
(def wd (fmat/decomposition-component truckdecomp :D))

(fmat/eigenvalues witch)
(fmat/eigenvalues truck)

(fmat/eigenvectors witch)
(fmat/eigenvectors truck)
(fmat/eigenvectors witch true)
(fmat/eigenvectors truck true)

;; ----------------------------------------


(def mdub (fm/seq->double-double-array [[1 2 3 4 5]
                                        [6 7 8 9 10]
                                        [-1 -2 -3 -4 -5]
                                        [-6 -7 -8 -9 -10]
                                        [1 2 3 4 5]]))
(type mdub)

(def m5x5 (fmat/rows->RealMatrix [[1 2 3 4 5]
                                  [6 7 8 9 10]
                                  [-1 -2 -3 -4 -5]
                                  [-6 -7 -8 -9 -10]
                                  [1 2 3 4 5]]))

;(fmat/eigenvectors m5x5)
; (err) Wrong number of args (5) passed to: fastmath.matrix/rows->mat

;(fmat/eigenvalues-matrix m5x5)
; (err) Wrong number of args (5) passed to: fastmath.matrix/rows->mat

(def m5x5-decomp (fmat/eigen-decomposition m5x5))

(fmat/decomposition-component m5x5-decomp :eigenvectors)
(fmat/decomposition-component (fmat/eigen-decomposition m5x5) :eigenvectors)
(map fvec/vec->seq
     (fmat/decomposition-component m5x5-decomp :eigenvectors))

;; These succeed:
(fmat/eigenvalues m5x5)
(fmat/decomposition-component m5x5-decomp :real-eigenvalues)
(fmat/decomposition-component m5x5-decomp :imag-eigenvalues)

(fmat/decomposition-component m5x5-decomp :V)
; #object[org.apache.commons.math3.linear.Array2DRowRealMatrix 0xc83fff9 "Array2DRowRealMatrix
; {{0.5379033159,-0.1179138184,1.4501010624,-1.0040216581,3.293351034},
;  {0.4882862549,-0.9871935558,-0.8123428453,1.7903179275,-7.7549758304},
;  {-0.5379033159,0.1179138184,-1.3828351958,-0.2704195191,3.4573201361},
;  {-0.4882862549,0.9871935558,-0.5977053223,-0.8140281116,3.1768830828},
;  {0.5379033159,-0.1179138184,1.3427823009,0.2981513614,-2.1725784226}}"]

(fmat/decomposition-component m5x5-decomp :D)
;#object[org.apache.commons.math3.linear.Array2DRowRealMatrix 0x276a9d93 "Array2DRowRealMatrix
;{{0.5,3.1224989992,0.0,0.0,0.0},
; {-3.1224989992,0.5,0.0,0.0,0.0},
; {0.0,0.0,0.0,0.0,0.0},
; {0.0,0.0,0.0,0.0,0.0},
; {0.0,0.0,0.0,0.0,0.0}}"]

