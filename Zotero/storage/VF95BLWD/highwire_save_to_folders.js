(function ($) {
  // Store our function as a property of Drupal.behaviors.
  Drupal.behaviors.highwire_saved_folder = { attach: function (context, settings) {
    var success = false;
    var id ='';
    var encodedUrl = '';
    
    if( $('.highwire_saved_search_link_ajax', context).length > 1){
      id = $('.highwire_saved_search_link_ajax').attr('id');
      if($('.highwire_saved_folder_form_ajax_'+id).length > 1) {
        $('.highwire_saved_folder_form_ajax_'+id+':first').remove();
      }
      $('select.ui-selectmenu-select').change(function(){
        $('div.form-item-folder-name').toggle($(this).val() == 'create_new_folder');
      });      
    }
    
    $('.highwire_saved_search_link_ajax', context).live( "click", function(){
      id = $(this).attr('id');
      $('.highwire_saved_folder_form_ajax_'+id).dialog({modal:true, draggable:false, width:540, title:'Save item in selected folder.'}); 
      return false;
    });
    
  $('input#edit-submit-add-to-folders').click(function(){
      if($('select#edit-folder-id').val() == 'create_new_folder')
      {
         if($('input#edit-folder-name').val() == '') {
            $('input#edit-folder-name').addClass('required error');
            $('input#edit-folder-name').parent().find('.description').html('Please enter a folder name.');   
            return false;      
         }
         if($('select#edit-folder-id').find(":contains('" + $.trim($('input#edit-folder-name').val()) + "')").length) {
            $('input#edit-folder-name').addClass('required error');
            $('input#edit-folder-name').parent().find('.description').html('Folder <strong>' + $('input#edit-folder-name').val()+ '</strong> already exists.' ); 
            return false;          
         } 
      }
      
      $('.ui-selectmenu-select').change(function(){
        if($(this).val() != 'create_new_folder') {
            $('input#edit-folder-name').removeClass('required error');
            $('input#edit-folder-name').parent().find('.description').html('Provide a name to new folder. e.g. <em>[keyword]</em>' );         
        }
      });
      
  }); 
    
  }};
}(jQuery));
