;; nvim functions and keymappings to evaluate scicloj Clay files

(module viz.lnvim
  {autoload {a aniseed.core
             str aniseed.string
             nvim aniseed.nvim
             eval conjure.eval
             extract conjure.extract}})

(var clay-not-yet-required true)

(defn require-clay-if-needed []
  "scicloj.clay.v2.api needs to be required in order for the eval-* functions
  to work, but this only needs to happen once per nrepl session.  This function
  ensures that this happens."
  (when clay-not-yet-required
    (set clay-not-yet-required false)
    (eval.eval-str 
      {:origin "custom-clay-wrapper"
       :code "(require 'scicloj.clay.v2.api)"})))

(defn eval-clojure-for-form-viz []
  "Display the result of the form that the cursor is in."
  (require-clay-if-needed)
  (eval.eval-str
    {:origin "custom-clay-wrapper"
     :code (str.join
             ""
             ["(scicloj.clay.v2.api/make! {:source-path \"" (nvim.fn.expand "%.") "\" :single-form `" (a.get (extract.form {:root? true}) :content) " :format [:html]})"])}))

;; In theory it would make sense to run the `require-clay-if-needed`
;; function below in the `on-filetype` function, and that does work.
;; However, it causes the the nrepl to be connected an extra time
;; and disconnected unnecessarily when on-filetype is first run.
;; Not sure why.  So instead I run `require-clay-if-needed` in the
;; keymapping callbacks.

(defn eval-clojure-for-ns-viz []
  "Display the entire namespace, typically a Clay source file."
  (require-clay-if-needed)
  (eval.eval-str
    {:origin "custom-clay-wrapper"
     :code (str.join
             ""
             ["(scicloj.clay.v2.api/make! {:source-path \"" (nvim.fn.expand "%.") "\" :format [:html]})"])}))

(defn on-filetype []
  "Map <localleader>ev to eval-clojure-for-form-viz and <localleader>env to
  eval-clojure-for-ns-viz."
  (nvim.buf_set_keymap
    0 :n "<localleader>ev" ""
    {:callback eval-clojure-for-form-viz})
  (nvim.buf_set_keymap
    0 :n "<localleader>env" ""
    {:callback eval-clojure-for-ns-viz}) )

(def augroup (nvim.create_augroup :viz.lnvim {}))
(nvim.create_autocmd
  :Filetype
  {:group augroup
   :pattern ["clojure"]
   :callback on-filetype})
