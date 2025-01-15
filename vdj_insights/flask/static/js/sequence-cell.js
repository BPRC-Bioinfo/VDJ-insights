document.addEventListener('DOMContentLoaded', function() {
    const sequenceCells = document.querySelectorAll('.sequence-cell');
    const modal = document.getElementById('sequenceModal');
    const modalHeader = document.getElementById('modalHeader');
    const modalContent = document.getElementById('fullSequence');
    const closeBtn = document.querySelector('.modal .close');
    const backdrop = document.querySelector('.modal-backdrop');

    sequenceCells.forEach(cell => {
        cell.addEventListener('click', function() {
            const sequence = this.getAttribute('data-sequence');
            const reference = this.closest('tr').querySelector('.reference-cell').textContent;
            modalHeader.textContent = reference; // Set the modal title to the Reference value
            modalContent.textContent = sequence;
            modal.style.display = 'block';
            backdrop.style.display = 'block';
        });
    });

    closeBtn.addEventListener('click', function() {
        modal.style.display = 'none';
        backdrop.style.display = 'none';
    });

    window.addEventListener('click', function(event) {
        if (event.target == backdrop) {
            modal.style.display = 'none';
            backdrop.style.display = 'none';
        }
    });
});
