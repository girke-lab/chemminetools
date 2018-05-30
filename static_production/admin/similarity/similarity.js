function clone(obj){
    if(obj == null || typeof(obj) != 'object')
        return obj;

    var temp = new obj.constructor(); // changed (twice)
    for(var key in obj)
        temp[key] = clone(obj[key]);

    return temp;
}

var count_text;	// hold text information created in HTML; avoid hard-coded text
// update the notice of "awaiting xx selected compounds"
function update_count() {
	jQuery("#selectionnotice").html(count_text);
	jQuery("#num_selection").text(2 - selection.length);
	if (selection.length == 2) {
		out_of_sight_left(jQuery("#selectionnotice"));
		// calc similarity here
		document.forms['compoundsel'].compound[0].value = compound_registry[selection[0].name.substring("check".length)];
		document.forms['compoundsel'].compound[1].value = compound_registry[selection[1].name.substring("check".length)];
		jQuery("#compoundsel").ajaxSubmit(ajax_similarity_option);
		return;
	}
	jQuery("#selectionnotice").html(count_text);
	jQuery("#num_selection").text(2 - selection.length);
}
// automatically arrange compounds in page
function auto_arrange() {
	pos = jQuery("#compoundarea").offset();
	jQuery(".compound").each(function() {
		jQuery(this).css('left', pos.left).css('top', pos.top);
		pos.left += 215;
		if (pos.left > jQuery(window).width() - 215) {
			pos.left = jQuery("#compoundarea").offset().left;
			pos.top += 310;
		}
	});
	jQuery("#compoundarea").css('height', pos.top + 310);
}
// process each compounds loaded by ajax using this function
function process_compound(cmp, animate) {
	if (animate == true)
		jQuery(cmp).show().fadeTo(1,0).delay(100).fadeTo(500,0.7).fadeTo(500,1).fadeTo(500,0.7).fadeTo(500,1).fadeTo(500,0.7).fadeTo(500,1)
	// make it draggable
	jQuery(cmp).show().draggable()
	// add mouse over and out highlights
	.mouseover(function() {
		jQuery(this).addClass("compoundhover");
	})
	.mouseout(function() {
		jQuery(this).removeClass("compoundhover");
	})
	// actions when the user checks a compound
	.children('.compoundcheck').click(function() {
		if (! this.checked) {
			// show the notice
			into_sight_left("#selectionnotice");
			// reset the similarity box
			jQuery("#startmcs").text('Load MCS');
			out_of_sight_left("#similarity");
			jQuery("#mcs-result").empty();
			jQuery("#mcssimilarities").hide();
			if (selection[0] == this) {selection = selection.slice(1);};
			if (selection[1] == this) {selection = selection.slice(0, 1);};
			update_count();
			jQuery(this.parentNode).removeClass('selectedcompound');
			return true;
		}
		if (selection.length == 2) {
			alert("You can only select two compounds to perform similarity calculation!");
			return false;
		} else {
			selection.push(this);
			update_count();
			jQuery(this.parentNode).addClass('selectedcompound');
		}
	});
	// actions when users close a compound box
	$(cmp).children("a.deletecompound").click(function() {
		if (confirm('Delete this compound?')) {
			jQuery(this).parent().hide();	
			var checkbox = jQuery(this).parent().children("input.compoundcheck").get(0);
			if (checkbox.checked)
				checkbox.click();
		}
		return false;
	});
}
var template = "";
var next_compound_id = 0;
var selection = new Array();

// general options for ajax
var ajax_option = {
	dataType : 'json',
	data: {ajax: '1'},
	error : function(req, status, err) {
		alert("Error: I was not able to submit the data to the server or the server cannot process your data.");
	},
	timeout : 50000,
	success : function(x) {
		alert("Server replies with " + x);
	}
};

// options for compound addition
var ajax_compound_upload_option = clone(ajax_option);
var compound_registry = {}; // for md5 to SMILES mapping
ajax_compound_upload_option.error = function(req, status, err) {
	alert("Error: I was not able to submit the data to the server or the SMILES or SDF you supplied is not valid.");
};
ajax_compound_upload_option.success = function(x) {
	all_status = x.status;
	failed = x.failed;
	if (failed != 0)
		alert(failed + 'compound(s) cannot be processed and is(are) ignored');
	for (var i = 0; i < all_status.length; i ++) {
		// add each compound returned to the workspace
		xx = all_status[i];
		compound_registry[xx.md5] = xx.smiles;
		if (xx['url'] == undefined)
			jQuery("#compoundarea").append(
				template
				.replace(/@url/, xx.img)
				.replace(/@id/, next_compound_id)
				.replace(/@md5/g, xx.md5)
				.replace(/@name/, xx.title)
			);
		else
			jQuery("#compoundarea").append(
				template
				.replace(/@url/, xx.img)
				.replace(/@id/, next_compound_id)
				.replace(/@md5/g, xx.md5)
				.replace(/@name/, '<a href="' + xx.url + '">' + xx.title + '</a>')
			);
		process_compound("#compound" + next_compound_id, all_status.length < 5);
		next_compound_id ++;
	}
	auto_arrange();
	// scroll to the bottom of the panel
	if (all_status.length)
		jQuery.scrollTo(jQuery("#compound" + (next_compound_id - 1)), 2000);
	// notify the user
	jQuery("<span>" + all_status.length + " loaded; </span>").appendTo(jQuery("#loadstatus")).hide().fadeIn('slow').animate({'opaque':1}, 3000, function() {jQuery(this).remove()});
};

// option for ajaxified MCS calculation
var ajax_mcs_option = clone(ajax_option);
ajax_mcs_option.error = function(req, status, err) {
	alert("Error: The server was not able to perform the requested calculation.");
};
ajax_mcs_option.success = function(x) {
	jQuery("#mcssize").text(x.m_size);
	jQuery("#mcsmin").text(x.sim_min);
	jQuery("#mcsmax").text(x.sim_max);
	jQuery("#mcstanimoto").text(x.sim_tanimoto);
	jQuery("#mcssimilarities").slideDown(500);
	jQuery("#compoundarea").append(
		template
		.replace(/@url/, x.img)
		.replace(/@id/, next_compound_id)
		.replace(/@md5/g, x.md5)
		.replace(/@name/, x.title)
	);
	compound_registry[x.md5] = x.smiles;
	jQuery("#mcs-result")
	.hide()
	.html(
		mcs_template
		.replace(/@url/, x.img)
	).slideDown(500);
	jQuery("#startmcs").text('Loaded');
	// let "add to workspace" work
	var last_loaded_mcs = next_compound_id;
	jQuery("#mcs-result #addmcsresult").unbind('click');
	jQuery("#mcs-result #addmcsresult").click(function(){
		jQuery(this).hide();
		var _offset = jQuery("#mcs-result").offset();
		jQuery("#compound" + last_loaded_mcs).css('left', _offset.left).css('top', _offset.top);
		process_compound("#compound" + last_loaded_mcs, false);
		jQuery("#compound" + last_loaded_mcs).animate({left: (_offset.left + 320) + 'px'}, 600);
		return false;
	});
	next_compound_id ++;
};

// option for ajaxified ap similarity calculation
var ajax_similarity_option = clone(ajax_option);
ajax_similarity_option.data['func-override'] = 'ap';
ajax_similarity_option.error = function(req, status, err) {
	alert("Error: The server was not able to calculate the similarity.");
};
ajax_similarity_option.success = function(x) {
	jQuery("#mcs-result").hide();
	jQuery("#mcssimilarities").hide();
	into_sight_left(jQuery("#similarityvalue").text(x.sim).parent());
	out_of_sight_left("#selectionnotice");
};

var similarity_left = -370; // where the similarity info box is placed
var mcs_template;
jQuery(document).ready(function() {
	similarity_left = jQuery("#similarity").offset().left;
	count_text = jQuery("#selectionnotice").html();
	template = jQuery("#compoundarea").find("script").html();
	mcs_template = jQuery("#mcs-result").find("script").html();
	jQuery("#compoundarea").empty().show();
	jQuery("#mcs-result").empty();
	jQuery("#mcssimilarities").hide();
	jQuery('#uploadform1, #uploadform2').ajaxForm(ajax_compound_upload_option);
	out_of_sight_top("#uploadarea", 10);
	//into_sight_top("#uploadarea");
	jQuery("#closeuploadbtn").click(function() {
		out_of_sight_top("#uploadarea");
		into_sight_top("#uploadareactrl");
		return false;
	});
	jQuery("a#startmcs").click(function() {
		if (selection.length != 2) {
			alert("You must select exactly 2 compounds to start MCS calculation");
			return false;
		}
		jQuery("#compoundsel").ajaxSubmit(ajax_mcs_option);
		return false;
	});
	jQuery(document).bind("ajaxStart", function() {
		jQuery("#ajaxstatus div").text("Working...");
		jQuery("#ajaxstatus div").show();
	});
	jQuery(document).bind("ajaxStop", function() {
		jQuery("#ajaxstatus div").text("Done")
			.queue(function() {$(this).show();$(this).dequeue();})
			.delay(2000)
			.queue(function() {$(this).hide();$(this).dequeue();});
	});
	jQuery("a#closesimilarity").click(function() {
		out_of_sight_left("#similarity");
		return false;
	});
	jQuery("a.resetworkspace").click(function() {
		if (confirm("Are you sure you want to clear all compounds in the workspace?"))
			jQuery("#compoundarea").empty();
		return false;
	});
	// load preloaded compounds is there is any
	if (preload) {
		ajax_compound_upload_option.success({status:preload, failed:[]});
	}
});

function out_of_sight_top(elem, duration) {
	if (duration == undefined) duration = 300;
	if (! jQuery(elem).hasClass('outofsight'))
		jQuery(elem).animate({top: '-=' + jQuery(elem).outerHeight() + 'px'}, duration).addClass('outofsight');
}

function into_sight_top(elem, duration) {
	if (duration == undefined) duration = 300;
	if (jQuery(elem).hasClass('outofsight'))
		jQuery(elem).show().animate({top: '+=' + jQuery(elem).outerHeight() + 'px'}, duration, 'swing').removeClass('outofsight');
}

function out_of_sight_left(elem, duration) {
	if (duration == undefined) duration = 300;
	if (! jQuery(elem).hasClass('outofsight'))
		jQuery(elem).animate({left: '-=' + jQuery(elem).outerWidth() + 'px'}, duration).addClass('outofsight');
}

function into_sight_left(elem, duration) {
	if (duration == undefined) duration = 300;
	if (jQuery(elem).hasClass('outofsight'))
		jQuery(elem).show().animate({left: '+=' + jQuery(elem).outerWidth() + 'px'}, duration).removeClass('outofsight');
}
