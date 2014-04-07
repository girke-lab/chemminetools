$(document).ready(function() {
	$("[data-toggle=popover]").popover({trigger: 'manual', html: true,});
	var currentElementId = false;
	$("[data-toggle=popover]").hover(function () {
		var cid = $(this).html();
		var thisCompound = $(this);
		currentElementId = thisCompound;
		$.ajax({
		    type: 'GET',
		    url: '/compounds/cid_lookup',
		    dataType: 'json',
		    success: function(data) {perm_id = data.id},
		    data: { cid: cid },
		    async: false
		});
		if(perm_id != 'ERROR'){
			var img = new Image();
			img.src = '/compounds/' + perm_id + '/png'
			thisCompound.attr('data-content', img.outerHTML);
			img.onload = function() {
				if(currentElementId == thisCompound){
					thisCompound.popover('fixTitle').popover('show');
				}
			}
		} else {
			$(this).attr('data-content', 'Error: cid not in workbench (or duplicated)')
			$(this).popover('fixTitle').popover('show');
		}
	}, function() {
		currentElementId = false;
		$(this).popover('hide');
	});

	gotoCidLink = function(){
		var cid = $(this).html();
		$.ajax({
		    type: 'GET',
		    url: '/compounds/cid_lookup',
		    dataType: 'json',
		    success: function(data) {perm_id = data.id},
		    data: { cid: cid },
		    async: false
		});
		if(perm_id != 'ERROR'){
			var cidLink = '/compounds/' + perm_id + '/';
			window.location.href = cidLink;
		}
	}

	$('.cidLink').click(gotoCidLink);
});
