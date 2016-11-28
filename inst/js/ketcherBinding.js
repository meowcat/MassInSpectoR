$(document).on("click", "button.ketcher", function(evt) {

  // evt.target is the button that was clicked
  var el = $(evt.target);
  // Raise an event to signal that the value changed
  el.trigger("change");
});


function getKetcher()
{
var frame = null;
if ('frames' in window && 'ketcherFrame' in window.frames)
frame = window.frames['ketcherFrame'];
else
return null;
if ('window' in frame)
return frame.window.ketcher;
}

function getSmiles()
{
var ketcher = getKetcher();
if (ketcher)
return(ketcher.getSmiles());
}

function getMolfile()
{
var ketcher = getKetcher();
if (ketcher)
return(ketcher.getMolfile());
}



function setMol(data)
{
var ketcher = getKetcher();
if (ketcher)
return(ketcher.setMolecule(data));
}


(function() {

  var ketcherInputBinding = new Shiny.InputBinding();
  $.extend(ketcherInputBinding, {
    find: function(scope) {
      return $(scope).find(".ketcher");
    },
    getValue: function(el) {
      mf = getSmiles();
      //alert(mf);
      return mf;
    },
    subscribe: function(el, callback) {
      $(el).on("change.ketcherInputBinding", function(e) {
        callback();
      });
    },
      // Receive messages from the server.
    // Messages sent by updateUrlInput() are received by this function.
    receiveMessage: function(el, data) {
      if (data.hasOwnProperty('value'))
      {
        setMol(data.value);
      }
      //$(el).trigger('change');
    }
    
  });
  Shiny.inputBindings.register(ketcherInputBinding);
})();



