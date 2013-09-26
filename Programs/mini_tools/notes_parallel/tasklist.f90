MODULE tasklist

  IMPLICIT NONE

  INTEGER :: max_task = 300
  INTEGER :: taskIndex

CONTAINS
  INTEGER FUNCTION get_next_task()

    ! Check if we are out of tasks
    IF (taskIndex == MAX_TASK) THEN
       get_next_task = -1
    ELSE
       taskIndex = taskIndex + 1
       get_next_task = taskIndex
    end if
    RETURN
  END FUNCTION get_next_task

  subroutine process_task(index)
    integer, intent(in) :: index

  end subroutine process_task

END MODULE tasklist


PROGRAM taskQueue

  USE taskList, ONLY: taskIndex, process_task, get_next_task
  INTEGER :: myIndex

  myIndex = get_next_task()
  DO WHILE (myindex /= -1)
     print *, myindex
     CALL process_task(myIndex)
     myIndex = get_next_task()
  ENDDO
END PROGRAM taskQueue