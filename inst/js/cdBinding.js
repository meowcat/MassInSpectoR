$(document).on("click", "button.chemdoodle_sketcher", function(evt) {

  // evt.target is the button that was clicked
  var el = $(evt.target);
  // Raise an event to signal that the value changed
  el.trigger("change");
});


(function() {

  var cdInputBinding = new Shiny.InputBinding();
  $.extend(cdInputBinding, {

    find: function(scope) {
      return $(scope).find(".chemdoodle_sketcher");
    },
    getValue: function(el) {
      var mol = sketcher.getMolecule();
      var jsonmol = new ChemDoodle.io.JSONInterpreter().molTo(mol);
      return jsonmol;
    },
    subscribe: function(el, callback) {
      $(el).on("change.cdInputBinding", function(e) {
        callback();
      });
    },

      // Receive messages from the server.
    // Messages sent by updateUrlInput() are received by this function.
    receiveMessage: function(el, data) {
      if (data.hasOwnProperty('value'))
      {
        var mol = ChemDoodle.readMOL(data.value);
        sketcher.loadMolecule(mol);
      }
      //$(el).trigger('change');
    }
    
  });
  Shiny.inputBindings.register(cdInputBinding);
})();



