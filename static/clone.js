$('#butt').click(function(){

var $orginal = $('form #container');
var $cloned = $orginal.clone();

//get original selects into a jq object
var $originalSelects = $orginal.find('select');
$cloned.find('select').each(function(index, item) {
     //set new select to value of old select
     $(item).val( $originalSelects.eq(index).val() );

});
//get original textareas into a jq object
var $originalTextareas = $orginal.find('textarea');

$cloned.find('textarea').each(function(index, item) {
//set new textareas to value of old textareas
$(item).val($originalTextareas.eq(index).val());
});
$cloned.appendTo('#clonedItem');
});
