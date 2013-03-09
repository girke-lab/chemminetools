var structures_loaded = false;
$(document).ready(function() {
	$("#loadstructures").click(function () {
		if (!structures_loaded) {
			var buf = $(this).text();
			$(".compounds .compound .image").each(function() {
			$(this).html('<img src="http://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid='+ $(this).attr('id').substring(4)  + '"/>');
		});
			structures_loaded = true;
		}
		if ($(this).get(0).checked)
			$(".compounds .compound .image").show();
		else
			$(".compounds .compound .image").hide();

	});
	$("#show-details").click(function() {
		$("#details").show();
		return false;
	});
	var n = $("td.distance").length;
	var cutoff = parseInt(parseFloat($($("td.distance").get(n-1)).text().substring('similarity: '.length)) * 100);
	if (cutoff == 100) cutoff = 99;
	$("#pugform").get(0).cutoff.value = cutoff;
	$("#pugform").submit(function() {
		cutoff = parseInt($(this).get(0).cutoff.value);
		if (!(cutoff <= 99 && cutoff >0)) {
			alert("You must type a cutoff value between 1 and 99.");
			return false;
		}
		if (cutoff < 10) $(this).get(0).cutoff.value = $(this).get(0).cutoff.value * 10;
		$(this).children("span").children("input[type=submit]").get(0).value = "Verifying. Please wait...";
		$(this).children("span").children("input[type=submit]").get(0).disabled = true;
	});
    $("#startpugsearch").click(function() {
        $("#pugform").modal();
        return false;
    });
});
