local _2afile_2a = ".lnvim.fnl"
local _2amodule_name_2a = "viz.lnvim"
local _2amodule_2a
do
  package.loaded[_2amodule_name_2a] = {}
  _2amodule_2a = package.loaded[_2amodule_name_2a]
end
local _2amodule_locals_2a
do
  _2amodule_2a["aniseed/locals"] = {}
  _2amodule_locals_2a = (_2amodule_2a)["aniseed/locals"]
end
local autoload = (require("nvim-local-fennel.aniseed.autoload")).autoload
local a, eval, extract, nvim, str = autoload("aniseed.core"), autoload("conjure.eval"), autoload("conjure.extract"), autoload("aniseed.nvim"), autoload("aniseed.string")
do end (_2amodule_locals_2a)["a"] = a
_2amodule_locals_2a["eval"] = eval
_2amodule_locals_2a["extract"] = extract
_2amodule_locals_2a["nvim"] = nvim
_2amodule_locals_2a["str"] = str
local clay_not_yet_required = true
local function require_clay_if_needed()
  if clay_not_yet_required then
    clay_not_yet_required = false
    return eval["eval-str"]({origin = "custom-clay-wrapper", code = "(require 'scicloj.clay.v2.api)"})
  else
    return nil
  end
end
_2amodule_2a["require-clay-if-needed"] = require_clay_if_needed
local function eval_clojure_for_form_viz()
  require_clay_if_needed()
  return eval["eval-str"]({origin = "custom-clay-wrapper", code = str.join("", {"(scicloj.clay.v2.api/make! {:source-path \"", nvim.fn.expand("%."), "\" :single-form `", a.get(extract.form({["root?"] = true}), "content"), " :format [:html]})"})})
end
_2amodule_2a["eval-clojure-for-form-viz"] = eval_clojure_for_form_viz
local function eval_clojure_for_ns_viz()
  require_clay_if_needed()
  return eval["eval-str"]({origin = "custom-clay-wrapper", code = str.join("", {"(scicloj.clay.v2.api/make! {:source-path \"", nvim.fn.expand("%."), "\" :format [:html]})"})})
end
_2amodule_2a["eval-clojure-for-ns-viz"] = eval_clojure_for_ns_viz
local function on_filetype()
  nvim.buf_set_keymap(0, "n", "<localleader>ev", "", {callback = eval_clojure_for_form_viz})
  return nvim.buf_set_keymap(0, "n", "<localleader>env", "", {callback = eval_clojure_for_ns_viz})
end
_2amodule_2a["on-filetype"] = on_filetype
local augroup = nvim.create_augroup("viz.lnvim", {})
do end (_2amodule_2a)["augroup"] = augroup
nvim.create_autocmd("Filetype", {group = augroup, pattern = {"clojure"}, callback = on_filetype})
return _2amodule_2a