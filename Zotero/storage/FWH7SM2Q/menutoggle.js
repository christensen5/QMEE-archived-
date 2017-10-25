function closeMenus(){
    $('.openMenu')
        .removeClass('openMenu')
        .slideUp(200)
    ;
}
function openMenu($menu){
    var open = $menu.hasClass('openMenu');
    closeMenus();
    if(!open){
        $menu
            .addClass('openMenu')
            .slideDown(200)
        ;
    }
}
    
$(document).ready(function () {	
    
	$('body').click(function(){
		closeMenus();
	});
	
    $('#toggleMore').click(function () {
		openMenu($('#moremenu'));
        return false;
    });
    $('#toggleMore2').click(function () {
        openMenu($('#moremenu2'));
        return false;
    });
    $('#toggleMore3').click(function () {   
        openMenu($('#moremenu3'));
        return false;
    });
    
    
    $('#toggleAccountMenu').click(function () {   
		openMenu($('#accountmenu'));
        return false;
    });
    
    $('#search_box input').click(function () {
    	closeMenus();
    });
    
    
    
    
});